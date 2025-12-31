#' Replace the third column with PRS
#'
#' Merges PRS values into a phenotype/covariate file and replaces column 3 with PRS.
#' Output is written to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @param dis_cov_file Character (path) OR list. Path to phenotype/covariate file
#'   (FID IID + >=1 columns). If a list is provided, it must contain `dis_cov_file`
#'   or `covar_file`.
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

  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)
  on_missing <- match.arg(on_missing)

  # ---- Allow passing previous-step lists ----
  original_cov_source <- NULL
  original_prs_source <- NULL

  if (is.list(dis_cov_file)) {
    if (!is.null(dis_cov_file$meta)) original_cov_source <- dis_cov_file$meta$step
    if (!is.null(dis_cov_file$dis_cov_file)) {
      dis_cov_file <- dis_cov_file$dis_cov_file
    } else if (!is.null(dis_cov_file$covar_file)) {
      dis_cov_file <- dis_cov_file$covar_file
    } else if (!is.null(dis_cov_file$output_file)) {
      dis_cov_file <- dis_cov_file$output_file
    } else {
      .stop2("Input list does not contain 'dis_cov_file', 'covar_file', or 'output_file'.")
    }
  }

  if (is.list(prs_file)) {
    if (!is.null(prs_file$meta)) original_prs_source <- prs_file$meta$step
    if (!is.null(prs_file$prs_file)) {
      prs_file <- prs_file$prs_file
    } else {
      .stop2("Input list does not contain 'prs_file' element.")
    }
  }

  # ---- Validate paths ----
  if (!is.character(dis_cov_file) || length(dis_cov_file) != 1 || !nzchar(dis_cov_file))
    .stop2("`dis_cov_file` must be a non-empty character scalar path.")
  if (!file.exists(dis_cov_file))
    .stop2("Covariate file not found: ", dis_cov_file)

  if (!is.character(prs_file) || length(prs_file) != 1 || !nzchar(prs_file))
    .stop2("`prs_file` must be a non-empty character scalar path.")
  if (!file.exists(prs_file))
    .stop2("PRS file not found: ", prs_file)

  # ---- Create tmp dir ----
  tmp_dir <- file.path(tempdir(), paste0("replace_cov_prs_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  message(.bar())
  message("GCIM-GWEIS-Z: Step 3 - Replace Exposure with PRS")
  message(.bar())
  message("Output directory: ", tmp_dir)

  # ---- Read covariate file robustly ----
  read_tab <- function(path, header) {
    utils::read.table(
      path,
      header = header,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    )
  }

  covar <- tryCatch(read_tab(dis_cov_file, header = TRUE),
                    error = function(e) NULL)

  if (is.null(covar) || ncol(covar) < 2 || !any(c("FID","IID","EID") %in% names(covar))) {
    # Retry: treat file as headerless
    covar <- tryCatch(
      read_tab(dis_cov_file, header = FALSE),
      error = function(e) .stop2("Failed to read covariate file: ", dis_cov_file, "\n", e$message)
    )

    if (nrow(covar) == 0) .stop2("Covariate file has 0 rows: ", dis_cov_file)
    if (ncol(covar) < 3) {
      .stop2("Covariate file must have at least 3 columns (FID IID + at least one covariate). ",
             "Found ", ncol(covar), " columns in: ", dis_cov_file)
    }

    # Assign required names by position
    names(covar)[1:2] <- c("FID", "IID")
    if (is.na(names(covar)[3]) || names(covar)[3] == "") names(covar)[3] <- "COVAR1"
  }

  # Normalize EID -> IID if needed
  if ("EID" %in% names(covar) && !"IID" %in% names(covar)) {
    names(covar)[names(covar) == "EID"] <- "IID"
  }

  if (!all(c("FID", "IID") %in% names(covar))) {
    .stop2(
      "Covariate file must contain 'FID' and 'IID' columns (or EID which will be renamed).\n",
      "Actual columns: ", paste(names(covar), collapse = ", "), "\n",
      "File: ", dis_cov_file
    )
  }

  if (ncol(covar) < 3) {
    .stop2("Covariate file must contain at least 3 columns (FID, IID, and column 3). File: ", dis_cov_file)
  }

  col3_name <- names(covar)[3]
  n_samples_original <- nrow(covar)

  message("Original covariate file:")
  message("  Samples: ", n_samples_original)
  message("  Columns: ", ncol(covar))
  message("  Column 3 to replace: '", col3_name, "'")

  # ---- Read PRS file ----
  prs <- tryCatch(
    utils::read.table(
      prs_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read PRS file: ", prs_file, "\n", e$message)
  )

  if (nrow(prs) == 0) .stop2("PRS file has 0 rows: ", prs_file)

  req <- c("FID", "IID", "PRS")
  miss <- setdiff(req, names(prs))
  if (length(miss) > 0) {
    .stop2(
      "PRS file missing required columns: ", paste(miss, collapse = ", "),
      "\nExpected: FID, IID, PRS\nActual: ", paste(names(prs), collapse = ", "),
      "\nFile: ", prs_file
    )
  }

  # ---- Check duplicates in PRS ----
  key <- paste(prs$FID, prs$IID, sep = "___")
  if (anyDuplicated(key)) {
    dups <- prs[duplicated(key) | duplicated(key, fromLast = TRUE), c("FID", "IID")]
    dup_file <- file.path(tmp_dir, "prs_duplicates.txt")
    utils::write.table(dups, dup_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    .stop2("PRS file has duplicated FID/IID pairs. Duplicates saved to: ", dup_file)
  }

  # ---- Match PRS to covariate file ----
  cov_key <- paste(covar$FID, covar$IID, sep = "___")
  prs_map <- prs$PRS
  names(prs_map) <- key
  prs_vec <- unname(prs_map[cov_key])

  n_missing <- sum(is.na(prs_vec))
  missing_file <- file.path(tmp_dir, "samples_missing_prs.txt")

  if (n_missing > 0) {
    missing_samples <- covar[is.na(prs_vec), c("FID", "IID"), drop = FALSE]
    utils::write.table(missing_samples, missing_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    if (on_missing == "stop") {
      .stop2(
        "Missing PRS values for ", n_missing, " samples (",
        sprintf("%.1f%%", 100 * n_missing / n_samples_original), ").\n",
        "Samples saved to: ", missing_file, "\n",
        "Fix: make sure FID/IID match between covariate and PRS; or use on_missing='drop'."
      )
    } else {
      warning("Dropping ", n_missing, " samples with missing PRS. See: ", missing_file, call. = FALSE)
      keep <- !is.na(prs_vec)
      covar <- covar[keep, , drop = FALSE]
      prs_vec <- prs_vec[keep]
    }
  }

  if (nrow(covar) == 0) .stop2("No samples remaining after PRS matching.")

  # ---- Replace column 3 with PRS ----
  covar[[3]] <- prs_vec
  names(covar)[3] <- "PRS"

  # ---- Save output ----
  out_file <- file.path(tmp_dir, out_name)
  utils::write.table(covar, out_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  message(.bar())
  message("PRS replacement completed successfully!")
  message("  Original column 3: '", col3_name, "'")
  message("  New column 3: 'PRS'")
  message("  Output samples: ", nrow(covar), " / ", n_samples_original)
  message("  Output file: ", out_file)
  message(.bar())

  res <- list(
    output_file = out_file,
    tmp_dir = tmp_dir,
    n_samples = nrow(covar),
    n_cols = ncol(covar),
    col3_replaced = col3_name,
    files = list(
      covar_prs = out_file,
      missing_prs = if (file.exists(missing_file)) missing_file else NA_character_
    ),
    meta = list(
      step = "PRS_REPLACE",
      replaced_col = col3_name,
      new_col_name = "PRS",
      on_missing = on_missing,
      n_samples_original = n_samples_original,
      n_samples_final = nrow(covar),
      n_samples_dropped = n_missing,
      source_covar = original_cov_source,
      source_prs = original_prs_source,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}