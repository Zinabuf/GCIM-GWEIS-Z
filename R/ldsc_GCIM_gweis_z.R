#' Run LDSC on GCIM-GWEIS interaction summary statistics
#'
#' This function munges GCIM-formatted interaction summary statistics,
#' runs LDSC heritability estimation, and extracts the LDSC intercept
#' for use in GCIM-GWEIS-Z correction. All outputs are saved to a
#' temporary directory.
#'
#' @param ldsc_path Character. Path to the LDSC directory (containing ldsc.py).
#' @param munge_path Character. Path to munge_sumstats.py script.
#' @param gcim_file Character. GCIM-formatted interaction file generated
#'   by \code{q_gweis()} or \code{b_gweis()}.
#' @param hm3_snplist Character. Path to HapMap3 SNP list (w_hm3.snplist).
#' @param ref_ld_chr Character. Path prefix to reference LD scores
#'   (for example, "eur_w_ld_chr/" for files like eur_w_ld_chr/1.l2.ldscore.gz).
#' @param chunksize Integer. Chunk size for munging (default is 100000).
#'
#' @return A list containing LDSC outputs stored in a temporary directory:
#' \describe{
#'   \item{munged_sumstats}{Path to munged .sumstats.gz file}
#'   \item{munge_log}{Path to munge log file}
#'   \item{h2_log}{Path to LDSC heritability log file}
#'   \item{intercept}{Numeric LDSC intercept value}
#'   \item{intercept_file}{Path to file containing intercept value}
#'   \item{intercept_se}{Standard error of intercept (if available)}
#'   \item{tmp_dir}{Temporary directory used by the function}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Validates all input paths (Python scripts, GCIM file, reference data)
#'   \item Creates a temporary directory for all outputs
#'   \item Munges summary statistics to LDSC format
#'   \item Runs LDSC heritability estimation
#'   \item Extracts intercept and standard error from the LDSC log
#'   \item Saves intercept to file for downstream GCIM-GWEIS-Z correction
#' }
#'
#' The LDSC intercept is used to correct for confounding in the GCIM-GWEIS-Z
#' framework. An intercept close to 1 suggests minimal inflation.
#'
#' @export
run_ldsc_gcim <- function(ldsc_path,
                          munge_path,
                          gcim_file,
                          hm3_snplist,
                          ref_ld_chr,
                          chunksize = 100000) {

  ## ---- Input validation ----
  if (!dir.exists(ldsc_path)) {
    stop("LDSC directory not found: ", ldsc_path)
  }

  ldsc_script <- file.path(ldsc_path, "ldsc.py")
  if (!file.exists(ldsc_script)) {
    stop("ldsc.py not found in: ", ldsc_path)
  }

  if (!file.exists(munge_path)) {
    stop("munge_sumstats.py not found at: ", munge_path)
  }

  if (!file.exists(gcim_file)) {
    stop("GCIM file not found: ", gcim_file)
  }

  if (!file.exists(hm3_snplist)) {
    stop("HapMap3 SNP list not found: ", hm3_snplist)
  }

  test_ld_file <- paste0(ref_ld_chr, "1.l2.ldscore.gz")
  if (!file.exists(test_ld_file)) {
    warning(
      "Could not find reference LD file: ", test_ld_file,
      "\nMake sure ref_ld_chr points to the correct prefix"
    )
  }

  if (!is.numeric(chunksize) || length(chunksize) != 1 ||
      is.na(chunksize) || chunksize < 1000) {
    stop("chunksize must be a positive integer >= 1000")
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("ldsc_gcim_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }
  message("Output directory: ", tmp_dir)

  out_prefix <- file.path(tmp_dir, "ldsc_gcim")

  ## ---- Helper: run python and capture logs ----
  run_python <- function(args, log_path) {
    # Capture both stdout and stderr into the same log file
    out <- tryCatch(
      system2("python", args = args, stdout = TRUE, stderr = TRUE),
      error = function(e) {
        writeLines(paste("ERROR:", e$message), con = log_path)
        return(structure(character(0), status = 1))
      }
    )
    writeLines(out, con = log_path)

    status <- attr(out, "status")
    if (is.null(status)) status <- 0
    return(status)
  }

  ## ---- Step 1: Munge summary statistics ----
  message("Step 1 of 3: Munging summary statistics")

  munge_log <- paste0(out_prefix, "_munge.log")
  munge_args <- c(
    shQuote(munge_path),
    "--sumstats", shQuote(gcim_file),
    "--chunksize", as.character(chunksize),
    "--merge-alleles", shQuote(hm3_snplist),
    "--out", shQuote(out_prefix)
  )

  munge_exit <- run_python(munge_args, munge_log)
  if (munge_exit != 0) {
    warning("Munge command exited with non-zero status. Check log: ", munge_log)
  }

  sumstats_file <- paste0(out_prefix, ".sumstats.gz")
  if (!file.exists(sumstats_file)) {
    stop(
      "Munged sumstats file not found: ", sumstats_file,
      "\nCheck munge log at: ", munge_log
    )
  }

  ## ---- Count SNPs in munged file (streaming) ----
  n_snps <- NA_integer_
  n_snps <- tryCatch(
    {
      con <- gzfile(sumstats_file, open = "rt")
      on.exit(close(con), add = TRUE)
      # Count lines without loading entire file
      n_lines <- 0L
      while (length(readLines(con, n = 10000L, warn = FALSE)) > 0) {
        n_lines <- n_lines + 10000L
      }
      # The loop overshoots by fixed chunk size; recount properly with a safer approach:
      close(con)

      con2 <- gzfile(sumstats_file, open = "rt")
      on.exit(close(con2), add = TRUE)
      n_lines2 <- 0L
      repeat {
        x <- readLines(con2, n = 50000L, warn = FALSE)
        if (length(x) == 0) break
        n_lines2 <- n_lines2 + length(x)
      }
      # subtract header line
      max(0L, n_lines2 - 1L)
    },
    error = function(e) NA_integer_
  )

  if (!is.na(n_snps)) {
    message("  SNPs in munged file: ", n_snps)
  }

  ## ---- Step 2: Run LDSC heritability ----
  message("Step 2 of 3: Running LDSC heritability estimation")

  h2_log <- paste0(out_prefix, ".log")
  ldsc_args <- c(
    shQuote(ldsc_script),
    "--h2", shQuote(sumstats_file),
    "--ref-ld-chr", shQuote(ref_ld_chr),
    "--w-ld-chr", shQuote(ref_ld_chr),
    "--out", shQuote(out_prefix)
  )

  ldsc_exit <- run_python(ldsc_args, h2_log)
  if (ldsc_exit != 0) {
    warning("LDSC command exited with non-zero status. Check log: ", h2_log)
  }

  if (!file.exists(h2_log)) {
    stop("LDSC log file not found: ", h2_log)
  }

  ## ---- Step 3: Extract intercept and SE ----
  message("Step 3 of 3: Extracting LDSC intercept")

  log_lines <- tryCatch(
    readLines(h2_log, warn = FALSE),
    error = function(e) {
      stop("Failed to read LDSC log file: ", h2_log, "\n", e$message)
    }
  )

  intercept_line <- grep("^Intercept:", log_lines, value = TRUE)
  if (length(intercept_line) == 0) {
    stop(
      "Intercept not found in LDSC log: ", h2_log,
      "\nLDSC may have failed. Check the log file for errors."
    )
  }

  # Typical line format: "Intercept: 1.0234 (0.0056)"
  intercept_parts <- strsplit(intercept_line[1], "\\s+")[[1]]
  intercept <- suppressWarnings(as.numeric(intercept_parts[2]))

  if (is.na(intercept)) {
    stop("Failed to parse intercept from line: ", intercept_line[1])
  }

  intercept_se <- NA_real_
  se_match <- regmatches(intercept_line[1], regexpr("\\(([0-9.]+)\\)", intercept_line[1]))
  if (length(se_match) > 0 && nzchar(se_match)) {
    intercept_se <- suppressWarnings(as.numeric(gsub("[()]", "", se_match)))
  }

  ## ---- Save intercept ----
  intercept_file <- paste0(out_prefix, "_intercept.txt")

  if (!is.na(intercept_se)) {
    writeLines(
      c(paste("Intercept:", intercept),
        paste("SE:", intercept_se)),
      con = intercept_file
    )
  } else {
    writeLines(as.character(intercept), con = intercept_file)
  }

  ## ---- Summary output ----
  message("LDSC analysis completed successfully")
  message("  Munged SNPs: ", ifelse(!is.na(n_snps), n_snps, "unknown"))
  message("  LDSC intercept: ", round(intercept, 4))
  if (!is.na(intercept_se)) {
    message("  Intercept SE: ", round(intercept_se, 4))
  }
  message("  Results saved to: ", tmp_dir)

  if (intercept > 1.1) {
    message("  NOTE: Intercept > 1.1 suggests potential confounding")
  } else if (intercept >= 0.9 && intercept <= 1.1) {
    message("  NOTE: Intercept close to 1 suggests minimal inflation")
  }

  ## ---- Return results ----
  result <- list(
    munged_sumstats = sumstats_file,
    munge_log = munge_log,
    h2_log = h2_log,
    intercept = intercept,
    intercept_file = intercept_file,
    intercept_se = if (!is.na(intercept_se)) intercept_se else NULL,
    tmp_dir = tmp_dir
  )

  invisible(result)
}