#' @importFrom utils read.table write.table head
#' @importFrom stats complete.cases
NULL

#' Run LDSC on PLINK GWEIS interaction summary statistics
#'
#' This function extracts interaction results (e.g., ADDxCOVAR1) from a PLINK2
#' .glm file, writes an LDSC-ready sumstats table (SNP, A1, A2, N, Z),
#' munges with munge_sumstats.py, runs ldsc.py --h2, and parses the LDSC intercept.
#' All outputs are saved to a temporary directory.
#'
#' @param ldsc_path Character. Path to LDSC directory containing ldsc.py.
#' @param munge_path Character. Path to munge_sumstats.py.
#' @param glm_file Character. Path to PLINK .glm.linear or .glm.logistic* file.
#' @param interaction_term Character. e.g., "ADDxCOVAR1".
#' @param hm3_snplist Character. Path to HapMap3 SNP list (w_hm3.snplist).
#' @param ref_ld_chr Character. Path prefix to reference LD scores (e.g. "eur_w_ld_chr/").
#' @param w_ld_chr Character. Path prefix to weights LD scores (often same as ref_ld_chr).
#' @param chunksize Integer. Chunk size for munging.
#' @param python Character. Python executable (default "python").
#'
#' @return A list containing LDSC outputs:
#' \describe{
#'   \item{ldsc_sumstats}{Path to intermediate LDSC-ready sumstats (TSV)}
#'   \item{munged_sumstats}{Path to munged .sumstats.gz}
#'   \item{munge_log}{Path to munge log file}
#'   \item{h2_log}{Path to LDSC heritability log file}
#'   \item{intercept}{Numeric LDSC intercept value}
#'   \item{intercept_se}{Standard error of intercept (if available)}
#'   \item{intercept_file}{Path to file containing intercept value}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{n_snps}{Number of SNPs written to LDSC-ready sumstats}
#' }
#'
#' @export
run_ldsc_gcim <- function(ldsc_path,
                          munge_path,
                          glm_file,
                          interaction_term,
                          hm3_snplist,
                          ref_ld_chr,
                          w_ld_chr = ref_ld_chr,
                          chunksize = 100000,
                          python = "python") {

  ## ---- Input validation ----
  if (!dir.exists(ldsc_path)) stop("LDSC directory not found: ", ldsc_path)
  ldsc_script <- file.path(ldsc_path, "ldsc.py")
  if (!file.exists(ldsc_script)) stop("ldsc.py not found in: ", ldsc_path)

  if (!file.exists(munge_path)) stop("munge_sumstats.py not found at: ", munge_path)
  if (!file.exists(glm_file)) stop("PLINK GLM file not found: ", glm_file)
  if (!file.exists(hm3_snplist)) stop("HapMap3 SNP list not found: ", hm3_snplist)

  test_ld_file <- paste0(ref_ld_chr, "1.l2.ldscore.gz")
  if (!file.exists(test_ld_file)) {
    warning(
      "Could not find reference LD file: ", test_ld_file,
      "\nMake sure ref_ld_chr points to the correct prefix"
    )
  }

  if (!is.character(interaction_term) || length(interaction_term) != 1 || !nzchar(interaction_term)) {
    stop("interaction_term must be a non-empty string (e.g., 'ADDxCOVAR1').")
  }

  if (!is.numeric(chunksize) || length(chunksize) != 1 || is.na(chunksize) || chunksize < 1000) {
    stop("chunksize must be a positive integer >= 1000")
  }

  ## ---- Create temporary directory ----
  tmp_dir <- file.path(tempdir(), paste0("ldsc_gcim_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", tmp_dir)

  out_prefix <- file.path(tmp_dir, "ldsc_gcim")

  ## ---- Helper: run python and write combined log ----
  run_python <- function(args, log_path) {
    out <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      error = function(e) {
        writeLines(paste("ERROR:", e$message), con = log_path)
        return(structure(character(0), status = 1))
      }
    )
    writeLines(out, con = log_path)
    status <- attr(out, "status")
    if (is.null(status)) status <- 0
    status
  }

  ## ---- Step 1: Read GLM and create LDSC-ready sumstats ----
  message("Step 1 of 3: Preparing LDSC-ready sumstats from PLINK GLM")

  glm <- tryCatch(
    read.table(glm_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE),
    error = function(e) stop("Failed to read GLM file: ", glm_file, "\n", e$message)
  )

  if (!"TEST" %in% names(glm)) stop("Column 'TEST' not found in GLM file: ", glm_file)

  sub <- glm[glm$TEST == interaction_term, , drop = FALSE]
  if (nrow(sub) == 0) {
    stop(
      "No rows found for interaction_term = ", interaction_term,
      "\nAvailable TEST values include: ",
      paste(head(unique(glm$TEST), 20), collapse = ", ")
    )
  }

  # Required columns for LDSC-ready table
  core_needed <- c("ID", "A1", "P", "OBS_CT")
  miss_core <- setdiff(core_needed, names(sub))
  if (length(miss_core) > 0) {
    stop(
      "Missing required columns in GLM subset: ", paste(miss_core, collapse = ", "),
      "\nColumns available: ", paste(names(sub), collapse = ", ")
    )
  }

  # Allele2 column (A2/REF/ALT/etc.)
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(sub)) {
      a2_col <- cand
      break
    }
  }
  if (is.null(a2_col)) {
    stop(
      "Could not find A2 column (tried A2/REF/ALT/ALT1). Columns: ",
      paste(names(sub), collapse = ", ")
    )
  }

  # Compute Z: prefer Z_STAT (logistic), else T_STAT (linear), else BETA/SE
  if ("Z_STAT" %in% names(sub)) {
    z <- suppressWarnings(as.numeric(sub$Z_STAT))
  } else if ("T_STAT" %in% names(sub)) {
    z <- suppressWarnings(as.numeric(sub$T_STAT))
  } else if (all(c("BETA", "SE") %in% names(sub))) {
    beta <- suppressWarnings(as.numeric(sub$BETA))
    se   <- suppressWarnings(as.numeric(sub$SE))
    z <- beta / se
  } else {
    stop(
      "Could not find Z_STAT or T_STAT, and cannot compute Z (need BETA and SE). Columns: ",
      paste(names(sub), collapse = ", ")
    )
  }

  n <- suppressWarnings(as.numeric(sub$OBS_CT))

  ldsc_df <- data.frame(
    SNP = sub$ID,
    A1  = sub$A1,
    A2  = sub[[a2_col]],
    N   = n,
    Z   = z,
    stringsAsFactors = FALSE
  )

  # Drop invalid
  keep <- complete.cases(ldsc_df) & is.finite(ldsc_df$Z) & ldsc_df$N > 0
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " SNPs with missing/invalid N or Z.")
    ldsc_df <- ldsc_df[keep, , drop = FALSE]
  }
  if (nrow(ldsc_df) == 0) stop("No valid SNPs remain after filtering invalid N/Z.")

  ldsc_sumstats <- file.path(tmp_dir, "ldsc_ready_sumstats.tsv")
  write.table(
    ldsc_df,
    file = ldsc_sumstats,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  n_snps <- nrow(ldsc_df)
  message("  SNPs written for LDSC munging: ", n_snps)

  ## ---- Step 2: Munge sumstats ----
  message("Step 2 of 3: Munging summary statistics")

  munge_log <- paste0(out_prefix, "_munge.log")
  munge_args <- c(
    munge_path,
    "--sumstats", ldsc_sumstats,
    "--snp", "SNP",
    "--a1", "A1",
    "--a2", "A2",
    "--N-col", "N",
    "--signed-sumstats", "Z,0",
    "--chunksize", as.character(chunksize),
    "--merge-alleles", hm3_snplist,
    "--out", out_prefix
  )

  munge_exit <- run_python(munge_args, munge_log)
  if (munge_exit != 0) warning("Munge exited with non-zero status. Check log: ", munge_log)

  munged_sumstats <- paste0(out_prefix, ".sumstats.gz")
  if (!file.exists(munged_sumstats)) {
    stop(
      "Munged sumstats file not found: ", munged_sumstats,
      "\nCheck munge log at: ", munge_log
    )
  }

  ## ---- Step 3: Run LDSC h2 ----
  message("Step 3 of 3: Running LDSC heritability estimation")

  h2_log <- paste0(out_prefix, ".log")
  ldsc_args <- c(
    ldsc_script,
    "--h2", munged_sumstats,
    "--ref-ld-chr", ref_ld_chr,
    "--w-ld-chr", w_ld_chr,
    "--out", out_prefix
  )

  ldsc_exit <- run_python(ldsc_args, h2_log)
  if (ldsc_exit != 0) warning("LDSC exited with non-zero status. Check log: ", h2_log)
  if (!file.exists(h2_log)) stop("LDSC log file not found: ", h2_log)

  ## ---- Parse intercept ----
  log_lines <- readLines(h2_log, warn = FALSE)
  intercept_line <- grep("^Intercept:", log_lines, value = TRUE)
  if (length(intercept_line) == 0) {
    stop("Intercept not found in LDSC log: ", h2_log, "\nCheck log for LDSC errors.")
  }

  # Typical: "Intercept: 1.0234 (0.0056)"
  parts <- strsplit(intercept_line[1], "\\s+")[[1]]
  intercept <- suppressWarnings(as.numeric(parts[2]))
  if (is.na(intercept)) stop("Failed to parse intercept from line: ", intercept_line[1])

  intercept_se <- NA_real_
  se_match <- regmatches(intercept_line[1], regexpr("\\(([0-9.]+)\\)", intercept_line[1]))
  if (length(se_match) > 0 && nzchar(se_match)) {
    intercept_se <- suppressWarnings(as.numeric(gsub("[()]", "", se_match)))
  }

  intercept_file <- paste0(out_prefix, "_intercept.txt")
  if (!is.na(intercept_se)) {
    writeLines(c(paste("Intercept:", intercept), paste("SE:", intercept_se)), con = intercept_file)
  } else {
    writeLines(as.character(intercept), con = intercept_file)
  }

  message("LDSC completed")
  message("  Intercept: ", round(intercept, 4))
  if (!is.na(intercept_se)) message("  Intercept SE: ", round(intercept_se, 4))

  invisible(list(
    ldsc_sumstats = ldsc_sumstats,
    munged_sumstats = munged_sumstats,
    munge_log = munge_log,
    h2_log = h2_log,
    intercept = intercept,
    intercept_se = if (!is.na(intercept_se)) intercept_se else NULL,
    intercept_file = intercept_file,
    tmp_dir = tmp_dir,
    n_snps = n_snps
  ))
}