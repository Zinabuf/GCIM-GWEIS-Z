#' LDSC intercept based Z score correction for GWEIS interaction results
#'
#' Applies LDSC intercept correction to interaction test statistics:
#'   Z_adj = Z_original / sqrt(intercept)
#' and recomputes two-sided p-values.
#'
#' @importFrom utils read.table write.table
#' @param glm_file Character. PLINK GWEIS output file (.glm.linear or .glm.logistic*).
#' @param intercept_file Character. File containing LDSC intercept value from run_ldsc_gcim().
#' @param int_term Character. Interaction term to extract (e.g., "ADDxCOVAR1").
#' @param out_file Character. Output filename (default "gcim_z_adjusted.txt").
#'
#' @return A list with adjusted results and output paths.
#' @export
gcim_z_adjust <- function(glm_file,
                          intercept_file,
                          int_term = "ADDxCOVAR1",
                          out_file = "gcim_z_adjusted.txt") {

  ## ---- Input validation ----
  if (!file.exists(glm_file)) stop("GWEIS GLM file not found: ", glm_file)
  if (!file.exists(intercept_file)) stop("LDSC intercept file not found: ", intercept_file)

  ## ---- Temp dir ----
  tmp_dir <- file.path(tempdir(), paste0("gcim_z_adjust_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", tmp_dir)

  ## ---- Step 1: Read GLM ----
  message("Step 1 of 4: Reading GWEIS results")
  gweis_res <- tryCatch(
    utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE,
               comment.char = "", check.names = FALSE),
    error = function(e) stop("Failed to read GWEIS file: ", glm_file, "\n", e$message)
  )
  message("  Total rows in GLM: ", nrow(gweis_res))

  if (!"TEST" %in% names(gweis_res)) stop("Column 'TEST' not found in GWEIS results")

  ## ---- Step 2: Extract interaction term ----
  message("Step 2 of 4: Extracting interaction term: ", int_term)
  available_tests <- unique(gweis_res$TEST)
  a <- gweis_res[gweis_res$TEST == int_term, , drop = FALSE]

  if (nrow(a) == 0) {
    stop("No interaction term found: ", int_term, "\nAvailable interaction terms: ",
         paste(grep("^ADDxCOVAR", available_tests, value = TRUE), collapse = ", "))
  }
  message("  Interaction SNPs extracted: ", nrow(a))

  ## ---- Step 3: Read intercept ----
  message("Step 3 of 4: Reading LDSC intercept")
  intercept_lines <- tryCatch(
    readLines(intercept_file, warn = FALSE),
    error = function(e) stop("Failed to read intercept file: ", intercept_file, "\n", e$message)
  )

  intercept_val <- suppressWarnings(as.numeric(trimws(intercept_lines[1])))

  # Fallback if file is like "Intercept: 1.02"
  if (is.na(intercept_val)) {
    idx <- grep("Intercept:", intercept_lines, ignore.case = TRUE)
    if (length(idx) > 0) {
      parts <- strsplit(intercept_lines[idx[1]], "\\s+")[[1]]
      if (length(parts) >= 2) intercept_val <- suppressWarnings(as.numeric(parts[2]))
    }
  }

  if (!is.finite(intercept_val) || intercept_val <= 0) {
    stop("Invalid LDSC intercept value in file: ", intercept_file)
  }
  message("  LDSC intercept: ", round(intercept_val, 4))

  ## ---- Step 4: Compute Z and adjust ----
  message("Step 4 of 4: Applying Z score correction")

  # Determine original Z statistic source
  z_source <- NULL
  if ("Z_STAT" %in% names(a)) {
    z <- suppressWarnings(as.numeric(a$Z_STAT))
    z_source <- "Z_STAT"
  } else if ("T_STAT" %in% names(a)) {
    z <- suppressWarnings(as.numeric(a$T_STAT))
    z_source <- "T_STAT"
  } else if (all(c("BETA", "SE") %in% names(a))) {
    beta <- suppressWarnings(as.numeric(a$BETA))
    se <- suppressWarnings(as.numeric(a$SE))
    z <- beta / se
    z_source <- "BETA/SE"
  } else {
    stop("Cannot compute Z: need Z_STAT or T_STAT or (BETA and SE). Columns: ",
         paste(names(a), collapse = ", "))
  }
  message("  Using Z source: ", z_source)

  # Original p-values
  if ("P" %in% names(a)) {
    p_orig <- suppressWarnings(as.numeric(a$P))
  } else {
    # If P absent, compute from Z
    p_orig <- 2 * stats::pnorm(-abs(z))
  }

  # Filter invalid
  keep <- is.finite(z) & !is.na(p_orig)
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " SNPs with missing/invalid Z or P.")
    a <- a[keep, , drop = FALSE]
    z <- z[keep]
    p_orig <- p_orig[keep]
  }
  if (nrow(a) == 0) stop("No valid SNPs remaining after filtering invalid statistics.")

  # Apply correction
  corr_factor <- sqrt(intercept_val)
  z_adj <- z / corr_factor
  p_adj <- 2 * stats::pnorm(-abs(z_adj))

  a$z_original <- z
  a$z_adj_int <- z_adj
  a$p_original <- p_orig
  a$p_value_int_adj <- p_adj

  n_sig_original <- sum(a$p_original < 0.05, na.rm = TRUE)
  n_sig_adjusted <- sum(a$p_value_int_adj < 0.05, na.rm = TRUE)

  message("  SNPs with p < 0.05 (original): ", n_sig_original)
  message("  SNPs with p < 0.05 (adjusted): ", n_sig_adjusted)

  ## ---- Save output ----
  output_path <- file.path(tmp_dir, out_file)
  utils::write.table(a, file = output_path, row.names = FALSE, quote = FALSE, sep = "\t")

  message("Z correction completed")
  message("  Total SNPs: ", nrow(a))
  message("  Intercept used: ", round(intercept_val, 4))
  message("  Correction factor sqrt(intercept): ", round(corr_factor, 4))
  message("  Saved: ", output_path)

  invisible(list(
    adjusted_data = a,
    output_file = output_path,
    tmp_dir = tmp_dir,
    intercept_used = intercept_val,
    correction_factor = corr_factor,
    z_source = z_source,
    n_snps = nrow(a),
    n_sig_original = n_sig_original,
    n_sig_adjusted = n_sig_adjusted
  ))
}