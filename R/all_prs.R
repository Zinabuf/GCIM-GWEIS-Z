#' Compute polygenic risk scores (PRS) in a target sample
#'
#' Computes PRS using PLINK2 based on SNP effect sizes (ID, A1, BETA)
#' and standardizes the resulting score (mean 0, SD 1) for downstream analyses.
#' All outputs are written to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @importFrom stats median sd
#' @param plink_path Character. Path to PLINK2 executable.
#' @param tar_mydata Character. Prefix of PLINK binary files (.bed/.bim/.fam).
#' @param score_file Character OR list. Either:
#'   (a) path to score file with columns: ID, A1, BETA; OR
#'   (b) a result list returned by \code{q_gwas()} or \code{b_gwas()} containing \code{$score_file}.
#' @param out_prefix Character. Prefix for PLINK output inside tmp_dir (default "prs").
#' @param threads Integer. Optional threads for PLINK (default 40).
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
#'   \item{files}{Named list of key output files (pipeline-friendly)}
#'   \item{meta}{Named list of metadata (pipeline-friendly)}
#' }
#'
#' @export
prs_scores <- function(plink_path,
                       tar_mydata,
                       score_file,
                       out_prefix = "prs",
                       threads = 40) {

  # ---- Helper functions ----
  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Resolve score_file if a previous-step list was provided ----
  gwas_source <- NULL
  if (is.list(score_file)) {
    if (!is.null(score_file$score_file)) {
      gwas_source <- if (!is.null(score_file$meta$trait_type)) {
        score_file$meta$trait_type
      } else {
        "unknown"
      }
      score_file <- score_file$score_file
    } else {
      .stop2("Input list does not contain a 'score_file' element.")
    }
  }

  # ---- Input validation ----
  if (!is.character(plink_path) || length(plink_path) != 1 || !nzchar(plink_path)) {
    .stop2("`plink_path` must be a non-empty character scalar.")
  }
  if (!file.exists(plink_path)) {
    .stop2("PLINK executable not found at: ", plink_path)
  }

  if (!is.character(tar_mydata) || length(tar_mydata) != 1 || !nzchar(tar_mydata)) {
    .stop2("`tar_mydata` must be a non-empty character scalar prefix to .bed/.bim/.fam.")
  }
  for (ext in c(".bed", ".bim", ".fam")) {
    fp <- paste0(tar_mydata, ext)
    if (!file.exists(fp)) .stop2("Missing PLINK file: ", fp)
  }

  if (!is.character(score_file) || length(score_file) != 1 || !nzchar(score_file)) {
    .stop2("`score_file` must be a non-empty character scalar path.")
  }
  if (!file.exists(score_file)) {
    .stop2("Score file not found: ", score_file)
  }

  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 1) {
    warning("Invalid `threads` value. Setting to 1.", call. = FALSE)
    threads <- 1
  }
  threads <- as.integer(threads)

  # ---- Validate score file header ----
  score_head <- tryCatch(
    utils::read.table(
      score_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      nrows = 5,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read score file: ", score_file, "\n", e$message)
  )
  
  req <- c("ID", "A1", "BETA")
  miss <- setdiff(req, names(score_head))
  if (length(miss) > 0) {
    .stop2(
      "Score file missing required columns: ", paste(miss, collapse = ", "),
      "\nExpected columns: ID, A1, BETA",
      "\nActual columns: ", paste(names(score_head), collapse = ", "),
      "\nFile: ", score_file
    )
  }

  # ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("prs_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 2 - PRS Construction (Target Sample)")
  message(.bar())
  message("Output directory: ", tmp_dir)
  if (!is.null(gwas_source)) {
    message("GWAS source: ", gwas_source, " trait")
  }
  message("Running PLINK2 PRS scoring...")

  # ---- Run PLINK score computation ----
  args <- c(
    "--bfile", tar_mydata,
    "--score", score_file, "1", "2", "3", "header",
    "--threads", as.character(threads),
    "--out", out_pref
  )

  exit_code <- suppressWarnings(
    system2(plink_path, args = args, stdout = stdout_file, stderr = stderr_file)
  )

  # Check for log file existence
  if (!file.exists(log_file)) {
    warning("Expected PLINK log not found: ", log_file, call. = FALSE)
  }

  if (!is.null(exit_code) && exit_code != 0) {
    .stop2(
      "PLINK exited with non-zero status (", exit_code, ").\n",
      "Check:\n  log: ", log_file, "\n  stderr: ", stderr_file, "\n  stdout: ", stdout_file
    )
  }

  # ---- Read PRS output ----
  sscore_file <- paste0(out_pref, ".sscore")
  if (!file.exists(sscore_file)) {
    .stop2("PRS output file not found: ", sscore_file, "\nCheck PLINK log at: ", log_file)
  }

  prs_raw <- tryCatch(
    utils::read.table(
      sscore_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read PRS output from ", sscore_file, "\n", e$message)
  )

  if (nrow(prs_raw) == 0) {
    .stop2("PRS output file has 0 rows: ", sscore_file)
  }

  if (!all(c("FID", "IID") %in% names(prs_raw))) {
    .stop2(
      "Columns 'FID' and 'IID' not found in PRS output.",
      "\nAvailable columns: ", paste(names(prs_raw), collapse = ", "),
      "\nFile: ", sscore_file
    )
  }

  # ---- Identify PRS sum column robustly ----
  score_sum_col <- NULL
  if ("SCORE1_SUM" %in% names(prs_raw)) {
    score_sum_col <- "SCORE1_SUM"
  } else {
    # Fallback: first column that matches SCORE*_SUM pattern
    cand <- grep("^SCORE\\d+_SUM$", names(prs_raw), value = TRUE)
    if (length(cand) > 0) score_sum_col <- cand[1]
  }
  
  if (is.null(score_sum_col)) {
    .stop2(
      "Could not find a PRS sum column (e.g., SCORE1_SUM) in: ", sscore_file,
      "\nAvailable columns: ", paste(names(prs_raw), collapse = ", ")
    )
  }

  message("Using PRS column: ", score_sum_col)

  # ---- Extract number of variants used (approx) ----
  n_variants <- NA_real_
  allele_ct_col <- NULL
  
  if ("ALLELE_CT" %in% names(prs_raw)) {
    allele_ct_col <- "ALLELE_CT"
    n_variants <- stats::median(prs_raw$ALLELE_CT, na.rm = TRUE)
  } else if ("NMISS_ALLELE_CT" %in% names(prs_raw)) {
    allele_ct_col <- "NMISS_ALLELE_CT"
    n_variants <- stats::median(prs_raw$NMISS_ALLELE_CT, na.rm = TRUE)
  }

  # ---- Build and scale PRS ----
  prs_df <- data.frame(
    FID = prs_raw$FID,
    IID = prs_raw$IID,
    PRS_raw = prs_raw[[score_sum_col]],
    stringsAsFactors = FALSE
  )

  # Convert to numeric and check for missing values
  prs_df$PRS_raw <- suppressWarnings(as.numeric(prs_df$PRS_raw))
  n_missing <- sum(is.na(prs_df$PRS_raw))
  
  if (n_missing > 0) {
    warning("Found ", n_missing, " samples with missing PRS values; removing them.", call. = FALSE)
    prs_df <- prs_df[!is.na(prs_df$PRS_raw), , drop = FALSE]
  }
  
  if (nrow(prs_df) == 0) {
    .stop2("No valid PRS values after removing missing data.")
  }

  # ---- Standardize PRS (mean=0, SD=1) ----
  mean_raw <- mean(prs_df$PRS_raw, na.rm = TRUE)
  sd_raw <- stats::sd(prs_df$PRS_raw)
  
  if (is.na(sd_raw) || isTRUE(all.equal(sd_raw, 0))) {
    warning(
      "PRS_raw has zero or NA variance; setting standardized PRS to 0 for all samples.",
      call. = FALSE
    )
    prs_df$PRS <- 0
  } else {
    prs_df$PRS <- (prs_df$PRS_raw - mean_raw) / sd_raw
  }

  # ---- Save scaled PRS file (FID IID PRS) ----
  prs_file <- file.path(tmp_dir, "prs_scaled.txt")
  utils::write.table(
    prs_df[, c("FID", "IID", "PRS"), drop = FALSE],
    file = prs_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  if (!file.exists(prs_file)) {
    .stop2("Failed to create scaled PRS file: ", prs_file)
  }

  message(.bar())
  message("PRS computation completed successfully!")
  message("  Samples with PRS: ", nrow(prs_df))
  if (!is.na(n_variants)) {
    message("  Variants used (approx): ", round(n_variants))
  }
  message("  PRS mean (raw): ", sprintf("%.4f", mean_raw))
  message("  PRS SD (raw):   ", sprintf("%.4f", sd_raw))
  message("  Output file: ", prs_file)
  message(.bar())

  # ---- Return results (pipeline-friendly structure) ----
  res <- list(
    # Backwards compatible fields
    sscore_file = sscore_file,
    prs_file = prs_file,
    prs_df = prs_df,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_samples = nrow(prs_df),
    n_variants = n_variants,

    # Pipeline-friendly additions
    files = list(
      sscore = sscore_file,
      prs_scaled = prs_file,
      plink_log = log_file,
      stdout = stdout_file,
      stderr = stderr_file
    ),
    meta = list(
      step = "PRS",
      gwas_source = gwas_source,
      score_sum_col = score_sum_col,
      allele_ct_col = allele_ct_col,
      threads = threads,
      n_samples_total = nrow(prs_raw),
      n_samples_valid = nrow(prs_df),
      n_samples_missing = n_missing,
      n_variants_approx = n_variants,
      prs_mean_raw = mean_raw,
      prs_sd_raw = sd_raw,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}