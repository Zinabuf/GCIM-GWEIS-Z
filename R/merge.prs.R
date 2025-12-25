#' Replace the third column with PRS
#'
#' Merges PRS values into a phenotype/covariate file and replaces column 3 with PRS.
#' Output is written to a temporary directory.
#'
#' @importFrom utils read.table write.table head
#' @param dis_cov_file Character. Path to phenotype/covariate file (FID IID + >=1 columns).
#' @param prs_file Character. Path to PRS file with columns FID, IID, PRS.
#' @param out_name Character. Output filename (default: "tar_cov_prs.txt").
#'
#' @return A list with:
#' \describe{
#'   \item{output_file}{Path to the output file with PRS replacing column 3}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{n_samples}{Number of samples in the output}
#'   \item{n_cols}{Number of columns in the output}
#'   \item{col3_replaced}{Name of the original third column that was replaced}
#' }
#'
#' @export
replace_covariate_with_prs <- function(dis_cov_file,
                                       prs_file,
                                       out_name = "tar_cov_prs.txt") {

  ## ---- Input validation ----
  if (!file.exists(dis_cov_file)) stop("Covariate file not found: ", dis_cov_file)
  if (!file.exists(prs_file)) stop("PRS file not found: ", prs_file)

  ## ---- Temp dir ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("replace_cov_prs_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", tmp_dir)

  ## ---- Read covariate file ----
  message("Step 1 of 4: Reading covariate file")
  covar <- tryCatch(
    utils::read.table(dis_cov_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read covariate file: ", dis_cov_file, "\n", e$message)
  )

  if (ncol(covar) < 3) stop("Covariate file must contain at least 3 columns (FID, IID, and a third column)")
  if (!all(c("FID", "IID") %in% names(covar))) stop("Covariate file must contain 'FID' and 'IID' columns")

  col3_name <- names(covar)[3]
  message("  Column 3 to be replaced: '", col3_name, "'")
  message("  Original samples: ", nrow(covar))

  ## ---- Read PRS file ----
  message("Step 2 of 4: Reading PRS file")
  prs <- tryCatch(
    utils::read.table(prs_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read PRS file: ", prs_file, "\n", e$message)
  )

  req <- c("FID", "IID", "PRS")
  miss <- setdiff(req, names(prs))
  if (length(miss) > 0) stop("PRS file missing required columns: ", paste(miss, collapse = ", "))

  # Ensure uniqueness in PRS mapping
  key <- paste(prs$FID, prs$IID, sep = "___")
  if (anyDuplicated(key)) {
    stop("PRS file has duplicated FID/IID pairs. Please deduplicate before merging.")
  }

  ## ---- Match PRS to covar (preserve row order) ----
  message("Step 3 of 4: Matching PRS to covariates (preserving order)")
  cov_key <- paste(covar$FID, covar$IID, sep = "___")
  prs_map <- prs$PRS
  names(prs_map) <- key

  prs_vec <- unname(prs_map[cov_key])

  n_missing <- sum(is.na(prs_vec))
  if (n_missing > 0) {
    warning("Found ", n_missing, " samples with missing PRS values")

    missing_samples <- covar[is.na(prs_vec), c("FID", "IID")]
    missing_file <- file.path(tmp_dir, "samples_missing_prs.txt")
    utils::write.table(missing_samples, file = missing_file,
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

    stop("Missing PRS values after matching for ", n_missing, " samples\n",
         "Samples with missing PRS saved to: ", missing_file, "\n",
         "Check that FID and IID match between covariate and PRS files")
  }

  ## ---- Replace column 3 with PRS ----
  message("Step 4 of 4: Replacing column 3 with PRS")
  covar[[3]] <- prs_vec

  ## ---- Save output ----
  out_file <- file.path(tmp_dir, out_name)
  utils::write.table(covar, file = out_file,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  message("Covariate replacement completed")
  message("  Final samples: ", nrow(covar))
  message("  Columns: ", ncol(covar))
  message("  Column 3 ('", col3_name, "') replaced with PRS")
  message("  Output saved to: ", out_file)

  invisible(list(
    output_file = out_file,
    tmp_dir = tmp_dir,
    n_samples = nrow(covar),
    n_cols = ncol(covar),
    col3_replaced = col3_name
  ))
}