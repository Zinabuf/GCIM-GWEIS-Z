#' Compute polygenic risk scores (PRS) in a target sample
#'
#' Computes PRS using PLINK2 based on SNP effect sizes (ID, A1, BETA)
#' and standardizes the resulting score (mean 0, SD 1) for downstream analyses.
#' All outputs are written to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @param plink_path Character. Path to PLINK2 executable.
#' @param tar_mydata Character. Prefix of PLINK binary files (.bed/.bim/.fam).
#' @param score_file Character. Path to score file with columns: ID, A1, BETA.
#' @param out_prefix Character. Prefix for PLINK output inside tmp_dir (default "prs").
#' @param threads Integer. Optional threads for PLINK (default 1).
#'
#' @return A list containing:
#' \describe{
#'   \item{sscore_file}{Path to PLINK .sscore file}
#'   \item{prs_file}{Path to scaled PRS file (FID IID PRS)}
#'   \item{prs_df}{Data frame with FID, IID, PRS_raw, PRS}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_samples}{Number of samples with PRS computed}
#'   \item{n_variants}{Approx number of variants used (median allele count)}
#' }
#'
#' @export
prs_scores <- function(plink_path,
                       tar_mydata,
                       score_file,
                       out_prefix = "prs",
                       threads = 1) {

  ## ---- Input validation ----
  if (!file.exists(plink_path)) {
    stop("PLINK executable not found at: ", plink_path)
  }

  for (ext in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(tar_mydata, ext))) {
      stop("Missing PLINK file: ", paste0(tar_mydata, ext))
    }
  }

  if (!file.exists(score_file)) {
    stop("Score file not found: ", score_file)
  }

  ## ---- Validate score file format (light check) ----
  score_head <- tryCatch(
    utils::read.table(score_file,
               header = TRUE,
               stringsAsFactors = FALSE,
               nrows = 5,
               check.names = FALSE),
    error = function(e) {
      stop("Failed to read score file: ", score_file, "\n", e$message)
    }
  )

  required_cols <- c("ID", "A1", "BETA")
  missing_cols <- setdiff(required_cols, names(score_head))
  if (length(missing_cols) > 0) {
    stop("Score file missing required columns: ",
         paste(missing_cols, collapse = ", "),
         "\nExpected columns: ID, A1, BETA")
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("prs_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")

  message("Output directory: ", tmp_dir)
  message("Step 1 of 3: Computing PRS with PLINK")

  ## ---- Run PLINK score computation ----
  args <- c(
    "--bfile", tar_mydata,
    "--score", score_file, "1", "2", "3", "header",
    "--threads", as.character(threads),
    "--out", out_pref
  )

  exit_code <- suppressWarnings(system2(plink_path, args = args))
  if (!is.null(exit_code) && exit_code != 0) {
    warning("PLINK exited with non-zero status (", exit_code, "). Check: ", log_file)
  }

  ## ---- Locate and read PRS output ----
  message("Step 2 of 3: Reading PRS output")

  sscore_file <- paste0(out_pref, ".sscore")
  if (!file.exists(sscore_file)) {
    stop("PRS output file not found: ", sscore_file,
         "\nCheck PLINK log at: ", log_file)
  }

  prs_raw <- tryCatch(
    utils::read.table(sscore_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      stop("Failed to read PRS output from ", sscore_file, "\n", e$message)
    }
  )

  if (!all(c("FID", "IID") %in% names(prs_raw))) {
    stop("Columns 'FID' and 'IID' not found in PRS output. Columns are: ",
         paste(names(prs_raw), collapse = ", "))
  }

  # Most common PLINK2 field for sum of scores
  if (!"SCORE1_SUM" %in% names(prs_raw)) {
    stop("Column 'SCORE1_SUM' not found in PRS output. Columns are: ",
         paste(names(prs_raw), collapse = ", "),
         "\nIf your PLINK output uses a different SCORE column, update the function.")
  }

  ## ---- Extract variants used ----
  n_variants <- NA_real_
  if ("ALLELE_CT" %in% names(prs_raw)) {
    n_variants <- stats::median(prs_raw$ALLELE_CT, na.rm = TRUE)
  } else if ("NMISS_ALLELE_CT" %in% names(prs_raw)) {
    n_variants <- stats::median(prs_raw$NMISS_ALLELE_CT, na.rm = TRUE)
  }

  ## ---- Build and scale PRS ----
  message("Step 3 of 3: Scaling PRS")

  prs_df <- data.frame(
    FID = prs_raw$FID,
    IID = prs_raw$IID,
    PRS_raw = prs_raw$SCORE1_SUM,
    stringsAsFactors = FALSE
  )

  # Drop missing PRS
  n_missing <- sum(is.na(prs_df$PRS_raw))
  if (n_missing > 0) {
    warning("Found ", n_missing, " samples with missing PRS values; removing them.")
    prs_df <- prs_df[!is.na(prs_df$PRS_raw), , drop = FALSE]
  }
  if (nrow(prs_df) == 0) stop("No valid PRS values after removing missing data.")

  # Handle zero variance safely
  if (stats::sd(prs_df$PRS_raw) == 0) {
    warning("PRS_raw has zero variance; setting standardized PRS to 0 for all samples.")
    prs_df$PRS <- 0
  } else {
    prs_df$PRS <- as.numeric(scale(prs_df$PRS_raw))
  }

  ## ---- Save scaled PRS file (FID IID PRS) ----
  prs_file <- file.path(tmp_dir, "prs_scaled.txt")
  utils::write.table(prs_df[, c("FID", "IID", "PRS")],
              file = prs_file,
              row.names = FALSE,
              col.names = TRUE,
              sep = "\t",
              quote = FALSE)

  message("PRS computation completed")
  message("  Samples with PRS: ", nrow(prs_df))
  if (!is.na(n_variants)) message("  Variants used (approx): ", round(n_variants))
  message("  Output PRS file: ", prs_file)

  invisible(list(
    sscore_file = sscore_file,
    prs_file = prs_file,
    prs_df = prs_df,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_samples = nrow(prs_df),
    n_variants = n_variants
  ))
}
