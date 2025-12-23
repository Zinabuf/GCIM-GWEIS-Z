#' Replace the third column with PRS
#'
#' This function merges PRS values into a phenotype and covariate file and
#' replaces the third column (the quantitative phenotype or interaction covariate)
#' with the computed PRS. This is used in the GCIM-GWEIS-Z pipeline to prepare
#' data for gene by environment interaction testing.
#'
#' @param qp_dis_cov Character. Path to the original phenotype and covariate file.
#'   Must contain columns FID and IID, and at least one additional column.
#' @param prs_file Character. Path to PRS file from \code{prs_scores()}.
#'   Must contain columns FID, IID, and PRS.
#' @param out_name Character. Name for output file (default: "qp_dis_cov_prs.txt").
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
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Validates input files exist and have required columns
#'   \item Reads phenotype and covariate file and PRS file
#'   \item Merges datasets by FID and IID
#'   \item Replaces column 3 with PRS values
#'   \item Saves output to a temporary directory
#' }
#'
#' The output file maintains the original column structure but with PRS
#' substituted into position 3, ready for analysis.
#'
#' @export
replace_covariate_with_prs <- function(qp_dis_cov,
                                       prs_file,
                                       out_name = "qp_dis_cov_prs.txt") {

  ## ---- Input validation ----
  if (!file.exists(qp_dis_cov)) {
    stop("Covariate file not found: ", qp_dis_cov)
  }

  if (!file.exists(prs_file)) {
    stop("PRS file not found: ", prs_file)
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("replace_cov_prs_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  message("Output directory: ", tmp_dir)

  ## ---- Read covariate file ----
  message("Step 1 of 4: Reading covariate file")

  covar <- tryCatch(
    {
      read.table(qp_dis_cov,
                 header = TRUE,
                 stringsAsFactors = FALSE,
                 check.names = FALSE)
    },
    error = function(e) {
      stop("Failed to read covariate file: ", qp_dis_cov, "\n", e$message)
    }
  )

  ## ---- Validate covariate file structure ----
  if (ncol(covar) < 3) {
    stop("Covariate file must contain at least 3 columns (FID, IID, and a third column)")
  }

  if (!all(c("FID", "IID") %in% names(covar))) {
    stop("Covariate file must contain 'FID' and 'IID' columns")
  }

  col3_name <- names(covar)[3]
  message("  Column 3 to be replaced: '", col3_name, "'")
  message("  Original samples: ", nrow(covar))

  ## ---- Read PRS file ----
  message("Step 2 of 4: Reading PRS file")

  prs <- tryCatch(
    {
      read.table(prs_file,
                 header = TRUE,
                 stringsAsFactors = FALSE,
                 check.names = FALSE)
    },
    error = function(e) {
      stop("Failed to read PRS file: ", prs_file, "\n", e$message)
    }
  )

  ## ---- Validate PRS file structure ----
  required_cols <- c("FID", "IID", "PRS")
  missing_cols <- setdiff(required_cols, names(prs))
  if (length(missing_cols) > 0) {
    stop("PRS file missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  message("  PRS samples: ", nrow(prs))

  ## ---- Merge covariate and PRS data ----
  message("Step 3 of 4: Merging PRS with covariates")

  covar_prs <- merge(
    covar,
    prs[, c("FID", "IID", "PRS")],
    by = c("FID", "IID"),
    all.x = TRUE
  )

  ## ---- Check for missing PRS values ----
  n_missing_prs <- sum(is.na(covar_prs$PRS))

  if (n_missing_prs > 0) {
    warning("Found ", n_missing_prs, " samples with missing PRS values")

    missing_samples <- covar_prs[is.na(covar_prs$PRS), c("FID", "IID")]

    missing_file <- file.path(tmp_dir, "samples_missing_prs.txt")
    write.table(
      missing_samples,
      file = missing_file,
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t",
      quote = FALSE
    )

    stop("Missing PRS values after merge for ", n_missing_prs, " samples\n",
         "Samples with missing PRS saved to: ", missing_file, "\n",
         "Check that FID and IID match between covariate and PRS files")
  }

  ## ---- Check sample count ----
  if (nrow(covar_prs) != nrow(covar)) {
    warning("Row count changed after merge: ",
            nrow(covar), " -> ", nrow(covar_prs))
  }

  ## ---- Replace column 3 with PRS ----
  message("Step 4 of 4: Replacing column 3 with PRS")

  original_col3 <- covar_prs[[3]]
  covar_prs[[3]] <- covar_prs$PRS
  covar_prs$PRS <- NULL

  ## ---- Restore original column order and names ----
  covar_prs <- covar_prs[, names(covar)]

  ## ---- Basic verification ----
  if (isTRUE(all.equal(covar_prs[[3]], original_col3))) {
    warning("Column 3 may be unchanged after replacement. Check PRS values.")
  }

  ## ---- Save output ----
  out_file <- file.path(tmp_dir, out_name)

  write.table(
    covar_prs,
    file = out_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  ## ---- Summary output ----
  message("Covariate replacement completed successfully")
  message("  Final samples: ", nrow(covar_prs))
  message("  Columns: ", ncol(covar_prs))
  message("  Column 3 ('", col3_name, "') replaced with PRS")
  message("  PRS statistics in column 3:")
  message("    Mean: ", round(mean(covar_prs[[3]], na.rm = TRUE), 4))
  message("    SD: ", round(sd(covar_prs[[3]], na.rm = TRUE), 4))
  message("    Range: [",
          round(min(covar_prs[[3]], na.rm = TRUE), 4), ", ",
          round(max(covar_prs[[3]], na.rm = TRUE), 4), "]")
  message("  Output saved to: ", out_file)

  invisible(list(
    output_file = out_file,
    tmp_dir = tmp_dir,
    n_samples = nrow(covar_prs),
    n_cols = ncol(covar_prs),
    col3_replaced = col3_name
  ))
}
