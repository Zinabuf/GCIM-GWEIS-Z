#' Replace the third column with PRS
#'
#' Merges PRS values into a phenotype/covariate file and replaces column 3 with PRS.
#' Output is written to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @param dis_cov_file Character (path) OR list. Path to phenotype/covariate file
#'   (FID IID + >=1 columns). If a list is provided, it must contain `dis_cov_file`
#'   or `covar_file` (you can adapt this later if needed).
#' @param prs_file Character (path) OR list from `prs_scores()` containing `$prs_file`.
#' @param out_name Character. Output filename (default: "tar_cov_prs.txt").
#' @param on_missing Character. What to do if PRS is missing for some samples:
#'   "stop" (default) or "drop".
#'
#' @return A list with:
#' \describe{
#'   \item{output_file}{Path to the output file with PRS replacing column 3}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{n_samples}{Number of samples in the output}
#'   \item{n_cols}{Number of columns in the output}
#'   \item{col3_replaced}{Name of the original third column that was replaced}
#'   \item{files}{Named list of key output files (pipeline-friendly)}
#'   \item{meta}{Named list of metadata (pipeline-friendly)}
#' }
#'
#' @export
replace_covariate_with_prs <- function(dis_cov_file,
                                       prs_file,
                                       out_name = "tar_cov_prs.txt",
                                       on_missing = c("stop", "drop")) {

  .stop2 <- function(...) stop(paste0(...), call. = FALSE)
  on_missing <- match.arg(on_missing)

  # ---- allow passing previous-step lists ----
  if (is.list(dis_cov_file)) {
    if (!is.null(dis_cov_file$dis_cov_file)) dis_cov_file <- dis_cov_file$dis_cov_file
    if (!is.null(dis_cov_file$covar_file))   dis_cov_file <- dis_cov_file$covar_file
  }
  if (is.list(prs_file) && !is.null(prs_file$prs_file)) {
    prs_file <- prs_file$prs_file
  }

  # ---- Input validation ----
  if (!is.character(dis_cov_file) || length(dis_cov_file) != 1 || !file.exists(dis_cov_file)) {
    .stop2("Covariate file not found: ", dis_cov_file)
  }
  if (!is.character(prs_file) || length(prs_file) != 1 || !file.exists(prs_file)) {
    .stop2("PRS file not found: ", prs_file)
  }

  # ---- Temp dir ----
  tmp_dir <- file.path(tempdir(), paste0("replace_cov_prs_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Read covariate file ----
  covar <- tryCatch(
    utils::read.table(dis_cov_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) .stop2("Failed to read covariate file: ", dis_cov_file, "\n", e$message)
  )

  if (ncol(covar) < 3) .stop2("Covariate file must contain at least 3 columns (FID, IID, and column 3).")
  if (!all(c("FID", "IID") %in% names(covar))) .stop2("Covariate file must contain 'FID' and 'IID' columns.")

  col3_name <- names(covar)[3]

  # ---- Read PRS file ----
  prs <- tryCatch(
    utils::read.table(prs_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) .stop2("Failed to read PRS file: ", prs_file, "\n", e$message)
  )

  req <- c("FID", "IID", "PRS")
  miss <- setdiff(req, names(prs))
  if (length(miss) > 0) .stop2("PRS file missing required columns: ", paste(miss, collapse = ", "))

  # Ensure uniqueness in PRS mapping
  key <- paste(prs$FID, prs$IID, sep = "___")
  if (anyDuplicated(key)) {
    .stop2("PRS file has duplicated FID/IID pairs. Please deduplicate before merging.")
  }

  # ---- Match PRS to covar (preserve row order) ----
  cov_key <- paste(covar$FID, covar$IID, sep = "___")
  prs_map <- prs$PRS
  names(prs_map) <- key
  prs_vec <- unname(prs_map[cov_key])

  n_missing <- sum(is.na(prs_vec))
  missing_file <- file.path(tmp_dir, "samples_missing_prs.txt")

  if (n_missing > 0) {
    missing_samples <- covar[is.na(prs_vec), c("FID", "IID"), drop = FALSE]
    utils::write.table(missing_samples, file = missing_file,
                       row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    if (on_missing == "stop") {
      .stop2(
        "Missing PRS values after matching for ", n_missing, " samples.\n",
        "Samples with missing PRS saved to: ", missing_file, "\n",
        "Check that FID/IID match between covariate and PRS files."
      )
    } else {
      warning("Dropping ", n_missing, " samples with missing PRS. See: ", missing_file, call. = FALSE)
      keep <- !is.na(prs_vec)
      covar <- covar[keep, , drop = FALSE]
      prs_vec <- prs_vec[keep]
    }
  }

  # ---- Replace column 3 with PRS ----
  covar[[3]] <- prs_vec

  # ---- Save output ----
  out_file <- file.path(tmp_dir, out_name)
  utils::write.table(covar, file = out_file,
                     row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  invisible(list(
    output_file = out_file,
    tmp_dir = tmp_dir,
    n_samples = nrow(covar),
    n_cols = ncol(covar),
    col3_replaced = col3_name,
    files = list(covar_prs = out_file, missing_prs = if (file.exists(missing_file)) missing_file else NA_character_),
    meta = list(step = "PRS_REPLACE_COL3", replaced_col = col3_name, on_missing = on_missing)
  ))
}