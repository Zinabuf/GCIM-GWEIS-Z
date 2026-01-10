#' Run LDSC heritability estimation (h2) and parse intercept
#'
#' @param ldsc_path Character. Path to LDSC directory containing `ldsc.py`,
#'   OR path directly to `ldsc.py`.
#' @param munged_sumstats Character OR list. Path to munged `.sumstats.gz`, or list
#'   from `munge_ldsc_gcim()` containing `$munged_sumstats` or `$files$munged`.
#' @param ref_ld_chr Character. Path prefix to reference LD scores (e.g. "eur_w_ld_chr/").
#' @param w_ld_chr Character. Path prefix to weights LD scores (default = ref_ld_chr).
#' @param python Character. Python executable (default "python").
#' @param tmp_dir Character. Optional output directory. If NULL, a new temp directory is created.
#' @param out_prefix Character. Output prefix within tmp_dir (default "ldsc_gcim").
#' @param show_log_tail_on_error Logical. Print log tail if LDSC fails.
#'
#' @return A list containing LDSC outputs and parsed intercept values.
#' @export
ldsc_h2_gcim <- function(ldsc_path,
                         munged_sumstats,
                         ref_ld_chr,
                         w_ld_chr = ref_ld_chr,
                         python = "python",
                         tmp_dir = NULL,
                         out_prefix = "ldsc_gcim",
                         show_log_tail_on_error = TRUE) {

  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)
  .tail_file <- function(path, n = 160) {
    if (!file.exists(path)) return(character(0))
    utils::tail(readLines(path, warn = FALSE), n)
  }

  # ---- Resolve munged_sumstats if list provided ----
  munge_source <- NULL
  if (is.list(munged_sumstats)) {
    if (!is.null(munged_sumstats$meta)) munge_source <- munged_sumstats$meta$step

    if (!is.null(munged_sumstats$munged_sumstats)) {
      munged_sumstats <- munged_sumstats$munged_sumstats
    } else if (!is.null(munged_sumstats$files$munged)) {
      munged_sumstats <- munged_sumstats$files$munged
    } else {
      .stop2("Input list does not contain 'munged_sumstats' or 'files$munged'.")
    }
  }

  # ---- Resolve ldsc.py path ----
  if (!is.character(ldsc_path) || length(ldsc_path) != 1 || !nzchar(ldsc_path))
    .stop2("`ldsc_path` must be a non-empty string.")

  ldsc_script <- ldsc_path
  if (dir.exists(ldsc_path)) {
    ldsc_script <- file.path(ldsc_path, "ldsc.py")
  }
  if (!file.exists(ldsc_script)) .stop2("ldsc.py not found at: ", ldsc_script)

  # ---- Validate files/paths ----
  if (!file.exists(munged_sumstats)) .stop2("Munged sumstats not found: ", munged_sumstats)
  if (!is.character(ref_ld_chr) || !nzchar(ref_ld_chr)) .stop2("ref_ld_chr must be non-empty.")
  if (!is.character(w_ld_chr) || !nzchar(w_ld_chr)) .stop2("w_ld_chr must be non-empty.")
  if (!is.character(python) || !nzchar(python)) .stop2("python must be non-empty.")

  # Basic LD-score existence checks (non-fatal warnings)
  test_ref <- paste0(ref_ld_chr, "1.l2.ldscore.gz")
  if (!file.exists(test_ref)) warning("Reference LD score not found: ", test_ref, call. = FALSE)
  test_w <- paste0(w_ld_chr, "1.l2.ldscore.gz")
  if (!file.exists(test_w)) warning("Weights LD score not found: ", test_w, call. = FALSE)

  # ---- Output dir ----
  if (is.null(tmp_dir)) {
    tmp_dir <- file.path(tempdir(), paste0("ldsc_h2_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_prefix_path <- file.path(tmp_dir, out_prefix)
  h2_log <- paste0(out_prefix_path, ".log")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 6 - LDSC Heritability Estimation")
  message(.bar())
  message("Output directory: ", tmp_dir)
  if (!is.null(munge_source)) message("Munged sumstats source: ", munge_source)

  # ---- Run LDSC (stdout+stderr -> log) ----
  ldsc_args <- c(
    ldsc_script,
    "--h2", munged_sumstats,
    "--ref-ld-chr", ref_ld_chr,
    "--w-ld-chr", w_ld_chr,
    "--out", out_prefix_path
  )

  message("Running LDSC: ", python, " ", paste(ldsc_args, collapse = " "))

  status <- system2(python, args = ldsc_args, stdout = h2_log, stderr = h2_log)
  if (is.null(status)) status <- 0

  if (status != 0) {
    if (isTRUE(show_log_tail_on_error) && file.exists(h2_log)) {
      cat("\n---- ldsc.py log (tail) ----\n")
      cat(paste(.tail_file(h2_log, 200), collapse = "\n"), "\n")
    }
    .stop2("LDSC exited with status ", status, ". See log: ", h2_log)
  }

  # ---- Parse intercept ----
  log_lines <- readLines(h2_log, warn = FALSE)
  intercept_line <- grep("^Intercept:", log_lines, value = TRUE)
  if (length(intercept_line) == 0) {
    .stop2("Intercept not found in LDSC log. See: ", h2_log)
  }

  line <- intercept_line[1]
  intercept <- suppressWarnings(as.numeric(sub("^Intercept:\\s*([0-9.eE+-]+).*$", "\\1", line)))
  if (is.na(intercept)) .stop2("Failed to parse intercept from: ", line)

  intercept_se <- NA_real_
  m <- regexec("\\(([0-9.eE+-]+)\\)", line)
  mm <- regmatches(line, m)[[1]]
  if (length(mm) >= 2) intercept_se <- suppressWarnings(as.numeric(mm[2]))

  message(.bar())
  message("LDSC completed successfully!")
  message("  Intercept: ", intercept, if (!is.na(intercept_se)) paste0(" (SE ", intercept_se, ")") else "")
  message("  Log: ", h2_log)
  message(.bar())

  invisible(list(
    h2_log = h2_log,
    intercept = intercept,
    intercept_se = intercept_se,
    tmp_dir = tmp_dir,
    files = list(h2_log = h2_log),
    meta = list(step = "LDSC_H2", munged_sumstats = munged_sumstats,
                ref_ld_chr = ref_ld_chr, w_ld_chr = w_ld_chr, timestamp = Sys.time())
  ))
}