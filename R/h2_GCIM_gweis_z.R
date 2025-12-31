#' Run LDSC heritability estimation (h2) and parse intercept
#'
#' Runs `ldsc.py --h2` on a munged `.sumstats.gz` file and parses the
#' LDSC intercept (and standard error if present) from the output log.
#' All outputs are saved to a temporary directory (or a user-provided directory).
#'
#' @param ldsc_path Character. Path to LDSC directory containing `ldsc.py`.
#' @param munged_sumstats Character OR list. Path to munged `.sumstats.gz`, or list
#'   from `munge_ldsc_gcim()` containing `$munged_sumstats` or `$files$munged`.
#' @param ref_ld_chr Character. Path prefix to reference LD scores (e.g. `"eur_w_ld_chr/"`).
#' @param w_ld_chr Character. Path prefix to weights LD scores (default = `ref_ld_chr`).
#' @param python Character. Python executable (default `"python"`).
#' @param tmp_dir Character. Optional output directory. If `NULL`, uses a new temp directory.
#' @param out_prefix Character. Optional output prefix within `tmp_dir` (default `"ldsc_gcim"`).
#'
#' @return A list containing LDSC outputs and parsed intercept values.
#' @export
ldsc_h2_gcim <- function(ldsc_path,
                         munged_sumstats,
                         ref_ld_chr,
                         w_ld_chr = ref_ld_chr,
                         python = "python",
                         tmp_dir = NULL,
                         out_prefix = "ldsc_gcim") {

  # ---- Helper functions ----
  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Resolve inputs if lists are provided ----
  munge_source <- NULL
  
  if (is.list(munged_sumstats)) {
    if (!is.null(munged_sumstats$meta)) {
      munge_source <- munged_sumstats$meta$step
    }
    
    if (!is.null(munged_sumstats$munged_sumstats)) {
      munged_sumstats <- munged_sumstats$munged_sumstats
    } else if (!is.null(munged_sumstats$files$munged)) {
      munged_sumstats <- munged_sumstats$files$munged
    } else {
      .stop2("Input list does not contain 'munged_sumstats' or 'files$munged' element.")
    }
  }

  # ---- Input validation ----
  if (!is.character(ldsc_path) || length(ldsc_path) != 1 || !nzchar(ldsc_path)) {
    .stop2("`ldsc_path` must be a non-empty character scalar.")
  }
  if (!dir.exists(ldsc_path)) {
    .stop2("LDSC directory not found: ", ldsc_path)
  }
  
  ldsc_script <- file.path(ldsc_path, "ldsc.py")
  if (!file.exists(ldsc_script)) {
    .stop2("ldsc.py not found in: ", ldsc_path)
  }

  if (!is.character(munged_sumstats) || length(munged_sumstats) != 1 || !nzchar(munged_sumstats)) {
    .stop2("`munged_sumstats` must be a non-empty character scalar path.")
  }
  if (!file.exists(munged_sumstats)) {
    .stop2("Munged sumstats file not found: ", munged_sumstats)
  }

  if (!is.character(ref_ld_chr) || length(ref_ld_chr) != 1 || !nzchar(ref_ld_chr)) {
    .stop2("`ref_ld_chr` must be a non-empty character scalar path prefix.")
  }

  if (!is.character(w_ld_chr) || length(w_ld_chr) != 1 || !nzchar(w_ld_chr)) {
    .stop2("`w_ld_chr` must be a non-empty character scalar path prefix.")
  }

  if (!is.character(python) || length(python) != 1 || !nzchar(python)) {
    .stop2("`python` must be a non-empty character scalar.")
  }

  # ---- Check for LD score files ----
  test_ld_file <- paste0(ref_ld_chr, "1.l2.ldscore.gz")
  if (!file.exists(test_ld_file)) {
    warning(
      "Reference LD file not found: ", test_ld_file, "\n",
      "  Make sure ref_ld_chr points to the correct prefix.\n",
      "  Expected format: 'path/to/eur_w_ld_chr/' (with trailing slash)",
      call. = FALSE
    )
  }

  test_w_file <- paste0(w_ld_chr, "1.l2.ldscore.gz")
  if (w_ld_chr != ref_ld_chr && !file.exists(test_w_file)) {
    warning(
      "Weight LD file not found: ", test_w_file, "\n",
      "  Make sure w_ld_chr points to the correct prefix.",
      call. = FALSE
    )
  }

  # ---- Create temporary directory ----
  if (is.null(tmp_dir)) {
    tmp_dir <- file.path(
      tempdir(),
      paste0("ldsc_h2_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  }
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_prefix_path <- file.path(tmp_dir, out_prefix)
  h2_log <- paste0(out_prefix_path, ".log")
  stdout_file <- paste0(out_prefix_path, ".stdout")
  stderr_file <- paste0(out_prefix_path, ".stderr")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 6 - LDSC Heritability Estimation")
  message(.bar())
  message("Output directory: ", tmp_dir)
  if (!is.null(munge_source)) {
    message("Munged sumstats source: ", munge_source)
  }
  message("Running LDSC heritability estimation...")

  # ---- Run LDSC h2 ----
  ldsc_args <- c(
    ldsc_script,
    "--h2", munged_sumstats,
    "--ref-ld-chr", ref_ld_chr,
    "--w-ld-chr", w_ld_chr,
    "--out", out_prefix_path
  )

  out <- tryCatch(
    system2(python, args = ldsc_args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      err_msg <- paste("ERROR:", e$message)
      writeLines(err_msg, con = h2_log)
      structure(err_msg, status = 1)
    }
  )

  # Write output to log files
  writeLines(out, con = h2_log)
  writeLines(out, con = stdout_file)
  
  # Separate stderr if available
  stderr_lines <- grep("(Error|ERROR|Warning|WARNING)", out, value = TRUE)
  if (length(stderr_lines) > 0) {
    writeLines(stderr_lines, con = stderr_file)
  }

  status <- attr(out, "status")
  if (is.null(status)) status <- 0

  if (status != 0) {
    .stop2(
      "LDSC exited with non-zero status (", status, ").",
      "\nCheck log: ", h2_log
    )
  }

  if (!file.exists(h2_log)) {
    .stop2("LDSC log file not found: ", h2_log)
  }

  # ---- Parse intercept from log ----
  log_lines <- readLines(h2_log, warn = FALSE)
  
  # Look for intercept line
  intercept_line <- grep("^Intercept:", log_lines, value = TRUE)
  if (length(intercept_line) == 0) {
    .stop2(
      "Intercept not found in LDSC log: ", h2_log,
      "\nCheck log for LDSC errors or warnings."
    )
  }

  # Expected format: "Intercept: 1.0123 (0.0045)"
  line <- intercept_line[1]
  message("Found intercept line: ", line)
  
  # Parse intercept value
  intercept <- suppressWarnings(
    as.numeric(sub("^Intercept:\\s*([0-9.eE+-]+).*$", "\\1", line))
  )
  
  if (is.na(intercept)) {
    .stop2(
      "Failed to parse intercept from line: ", line,
      "\nExpected format: 'Intercept: <value> (<se>)'"
    )
  }

  # Parse standard error (if present)
  intercept_se <- NA_real_
  m <- regexec("\\(([0-9.eE+-]+)\\)", line)
  mm <- regmatches(line, m)[[1]]
  if (length(mm) >= 2) {
    intercept_se <- suppressWarnings(as.numeric(mm[2]))
  }

  # ---- Parse additional h2 statistics (optional) ----
  h2_line <- grep("^Total Observed scale h2:", log_lines, value = TRUE)
  h2_value <- NA_real_
  h2_se <- NA_real_
  
  if (length(h2_line) > 0) {
    h2_match <- regexec("([0-9.eE+-]+)\\s*\\(([0-9.eE+-]+)\\)", h2_line[1])
    h2_parts <- regmatches(h2_line[1], h2_match)[[1]]
    if (length(h2_parts) >= 3) {
      h2_value <- suppressWarnings(as.numeric(h2_parts[2]))
      h2_se <- suppressWarnings(as.numeric(h2_parts[3]))
    }
  }

  # ---- Write intercept files ----
  intercept_file <- paste0(out_prefix_path, "_intercept.txt")
  if (!is.na(intercept_se)) {
    writeLines(
      c(paste("Intercept:", intercept), paste("SE:", intercept_se)),
      con = intercept_file
    )
  } else {
    writeLines(as.character(intercept), con = intercept_file)
  }

  if (!file.exists(intercept_file)) {
    .stop2("Failed to create intercept file: ", intercept_file)
  }

  # ---- Write intercept TSV ----
  intercept_tsv <- paste0(out_prefix_path, "_intercept.tsv")
  intercept_df <- data.frame(
    intercept = intercept,
    intercept_se = if (!is.na(intercept_se)) intercept_se else NA_real_,
    stringsAsFactors = FALSE
  )
  
  utils::write.table(
    intercept_df,
    file = intercept_tsv,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  if (!file.exists(intercept_tsv)) {
    .stop2("Failed to create intercept TSV file: ", intercept_tsv)
  }

  message(.bar())
  message("LDSC heritability estimation completed successfully!")
  message("  Intercept: ", sprintf("%.6f", intercept))
  if (!is.na(intercept_se)) {
    message("  Intercept SE: ", sprintf("%.6f", intercept_se))
  }
  if (!is.na(h2_value)) {
    message("  h2 (observed scale): ", sprintf("%.6f", h2_value))
    if (!is.na(h2_se)) {
      message("  h2 SE: ", sprintf("%.6f", h2_se))
    }
  }
  message("  Log file: ", h2_log)
  message("  Intercept file: ", intercept_file)
  message("  Intercept TSV: ", intercept_tsv)
  message(.bar())

  # ---- Return results (pipeline-friendly structure) ----
  res <- list(
    # Backwards compatible fields
    h2_log = h2_log,
    intercept = intercept,
    intercept_se = if (!is.na(intercept_se)) intercept_se else NULL,
    intercept_file = intercept_file,
    intercept_tsv = intercept_tsv,
    tmp_dir = tmp_dir,

    # Pipeline-friendly additions
    files = list(
      h2_log = h2_log,
      intercept_txt = intercept_file,
      intercept_tsv = intercept_tsv,
      stdout = stdout_file,
      stderr = if (file.exists(stderr_file)) stderr_file else NA_character_
    ),
    meta = list(
      step = "LDSC_H2",
      munged_sumstats = munged_sumstats,
      ref_ld_chr = ref_ld_chr,
      w_ld_chr = w_ld_chr,
      munge_source = munge_source,
      intercept = intercept,
      intercept_se = if (!is.na(intercept_se)) intercept_se else NA_real_,
      h2 = if (!is.na(h2_value)) h2_value else NA_real_,
      h2_se = if (!is.na(h2_se)) h2_se else NA_real_,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}