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
#' @return A list containing LDSC outputs.
#' @export
ldsc_h2_gcim <- function(ldsc_path,
                         munged_sumstats,
                         ref_ld_chr,
                         w_ld_chr = ref_ld_chr,
                         python = "python",
                         tmp_dir = NULL,
                         out_prefix = "ldsc_gcim") {

  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # allow passing prior-step lists
  if (is.list(munged_sumstats) && !is.null(munged_sumstats$munged_sumstats)) {
    munged_sumstats <- munged_sumstats$munged_sumstats
  }
  if (is.list(munged_sumstats) && !is.null(munged_sumstats$files$munged)) {
    munged_sumstats <- munged_sumstats$files$munged
  }

  if (!dir.exists(ldsc_path)) .stop2("LDSC directory not found: ", ldsc_path)
  ldsc_script <- file.path(ldsc_path, "ldsc.py")
  if (!file.exists(ldsc_script)) .stop2("ldsc.py not found in: ", ldsc_path)

  if (!is.character(munged_sumstats) || length(munged_sumstats) != 1 || !file.exists(munged_sumstats)) {
    .stop2("Munged sumstats file not found: ", munged_sumstats)
  }

  # light sanity check for LD score prefix
  test_ld_file <- paste0(ref_ld_chr, "1.l2.ldscore.gz")
  if (!file.exists(test_ld_file)) {
    warning("Could not find reference LD file: ", test_ld_file,
            "\nMake sure ref_ld_chr points to the correct prefix", call. = FALSE)
  }

  if (is.null(tmp_dir)) {
    tmp_dir <- file.path(tempdir(), paste0("ldsc_h2_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_prefix_path <- file.path(tmp_dir, out_prefix)
  h2_log <- paste0(out_prefix_path, ".log")
  stdout_file <- paste0(out_prefix_path, ".stdout")
  stderr_file <- paste0(out_prefix_path, ".stderr")

  message("Running LDSC heritability estimation (h2)")
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
      writeLines(paste("ERROR:", e$message), con = h2_log)
      structure(character(0), status = 1)
    }
  )
  writeLines(out, con = h2_log)

  status <- attr(out, "status"); if (is.null(status)) status <- 0
  # also keep copies
  writeLines(out, con = stdout_file)
  writeLines(out, con = stderr_file)

  if (status != 0) {
    .stop2("LDSC exited with non-zero status (", status, "). Check log: ", h2_log)
  }
  if (!file.exists(h2_log)) .stop2("LDSC log file not found: ", h2_log)

  # ---- Parse intercept ----
  log_lines <- readLines(h2_log, warn = FALSE)
  intercept_line <- grep("^Intercept:", log_lines, value = TRUE)
  if (length(intercept_line) == 0) {
    .stop2("Intercept not found in LDSC log: ", h2_log, "\nCheck log for LDSC errors.")
  }

  # Expected like: "Intercept: 1.0123 (0.0045)"
  line <- intercept_line[1]
  intercept <- suppressWarnings(as.numeric(sub("^Intercept:\\s*([0-9.eE+-]+).*$", "\\1", line)))
  if (is.na(intercept)) .stop2("Failed to parse intercept from line: ", line)

  intercept_se <- NA_real_
  m <- regexec("\\(([0-9.eE+-]+)\\)", line)
  mm <- regmatches(line, m)[[1]]
  if (length(mm) >= 2) {
    intercept_se <- suppressWarnings(as.numeric(mm[2]))
  }

  # ---- Write intercept files ----
  intercept_file <- paste0(out_prefix_path, "_intercept.txt")
  if (!is.na(intercept_se)) {
    writeLines(c(paste("Intercept:", intercept), paste("SE:", intercept_se)), con = intercept_file)
  } else {
    writeLines(as.character(intercept), con = intercept_file)
  }

  intercept_tsv <- paste0(out_prefix_path, "_intercept.tsv")
  utils::write.table(
    data.frame(intercept = intercept, intercept_se = if (!is.na(intercept_se)) intercept_se else NA_real_),
    file = intercept_tsv, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
  )

  invisible(list(
    h2_log = h2_log,
    intercept = intercept,
    intercept_se = if (!is.na(intercept_se)) intercept_se else NULL,
    intercept_file = intercept_file,
    intercept_tsv = intercept_tsv,
    tmp_dir = tmp_dir,
    files = list(
      h2_log = h2_log,
      intercept_txt = intercept_file,
      intercept_tsv = intercept_tsv
    ),
    meta = list(
      step = "LDSC_H2",
      munged_sumstats = munged_sumstats,
      ref_ld_chr = ref_ld_chr,
      w_ld_chr = w_ld_chr
    )
  ))
}