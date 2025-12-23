#' LDSC intercept based Z score correction for GWEIS interaction results
#'
#' This function applies LDSC intercept correction to SNP by covariate
#' interaction test statistics and computes adjusted Z scores and p values.
#' All outputs are saved to a temporary directory.
#'
#' @param glm_file Character. PLINK GWEIS output file
#'   (.glm.linear or .glm.logistic.hybrid) from \code{q_gweis()} or \code{b_gweis()}.
#' @param intercept_file Character. File containing the LDSC intercept value
#'   from \code{run_ldsc_gcim()}.
#' @param int_term Character. Interaction term name to extract
#'   (default is "ADDxCOVAR1").
#' @param out_file Character. Output file name
#'   (default is "gcim_z_adjusted.txt").
#'
#' @return A list containing adjusted results stored in a temporary directory:
#' \describe{
#'   \item{adjusted_data}{Data frame with adjusted Z scores and p values}
#'   \item{output_file}{Path to saved adjusted results file}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{intercept_used}{LDSC intercept value used for correction}
#'   \item{correction_factor}{Square root of the intercept used as the divisor}
#'   \item{n_snps}{Number of SNPs in adjusted results}
#'   \item{n_sig_original}{Number of significant SNPs before correction (p < 0.05)}
#'   \item{n_sig_adjusted}{Number of significant SNPs after correction (p < 0.05)}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Reads PLINK GWEIS results and extracts the interaction term
#'   \item Reads the LDSC intercept from file
#'   \item Applies Z score correction: Z_adj = Z_original / sqrt(intercept)
#'   \item Computes two sided adjusted p values
#'   \item Saves corrected results to a temporary directory
#' }
#'
#' The correction adjusts for inflation using the LDSC intercept. The adjusted
#' Z score is computed as Z_adj = Z_original / sqrt(intercept).
#'
#' @export
gcim_z_adjust <- function(glm_file,
                          intercept_file,
                          int_term = "ADDxCOVAR1",
                          out_file = "gcim_z_adjusted.txt") {

  ## ---- Input validation ----
  if (!file.exists(glm_file)) {
    stop("GWEIS GLM file not found: ", glm_file)
  }

  if (!file.exists(intercept_file)) {
    stop("LDSC intercept file not found: ", intercept_file)
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("gcim_z_adjust_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }
  message("Output directory: ", tmp_dir)

  ## ---- Step 1: Read GWEIS results ----
  message("Step 1 of 4: Reading GWEIS results")

  gweis_res <- tryCatch(
    {
      read.table(
        glm_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        comment.char = "",
        check.names = FALSE
      )
    },
    error = function(e) {
      stop("Failed to read GWEIS file: ", glm_file, "\n", e$message)
    }
  )

  message("  Total rows in GWEIS file: ", nrow(gweis_res))

  if (!"TEST" %in% colnames(gweis_res)) {
    stop("Column 'TEST' not found in GWEIS results")
  }

  required_cols <- c("T_STAT", "P")
  missing_cols <- setdiff(required_cols, colnames(gweis_res))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in GWEIS results: ",
         paste(missing_cols, collapse = ", "))
  }

  ## ---- Step 2: Extract interaction term ----
  message("Step 2 of 4: Extracting interaction term: ", int_term)

  available_tests <- unique(gweis_res$TEST)
  message("  Available TEST values: ", paste(available_tests, collapse = ", "))

  a <- gweis_res[gweis_res$TEST == int_term, ]

  if (nrow(a) == 0) {
    stop(
      "No interaction term found: ", int_term, "\n",
      "Available interaction terms: ",
      paste(grep("ADDx", available_tests, value = TRUE), collapse = ", ")
    )
  }

  message("  Interaction SNPs extracted: ", nrow(a))

  ## ---- Step 3: Read and validate LDSC intercept ----
  message("Step 3 of 4: Reading LDSC intercept")

  intercept_lines <- tryCatch(
    readLines(intercept_file, warn = FALSE),
    error = function(e) {
      stop("Failed to read LDSC intercept file: ", intercept_file, "\n", e$message)
    }
  )

  intercept_val <- NA_real_
  intercept_val <- suppressWarnings(as.numeric(trimws(intercept_lines[1])))

  if (is.na(intercept_val)) {
    idx <- grep("Intercept:", intercept_lines, ignore.case = TRUE)
    if (length(idx) > 0) {
      parts <- strsplit(intercept_lines[idx[1]], "\\s+")[[1]]
      if (length(parts) >= 2) {
        intercept_val <- suppressWarnings(as.numeric(parts[2]))
      }
    }
  }

  if (is.na(intercept_val) || !is.finite(intercept_val)) {
    stop("Invalid LDSC intercept value in file: ", intercept_file)
  }

  if (intercept_val <= 0) {
    stop("LDSC intercept must be positive, got: ", intercept_val)
  }

  message("  LDSC intercept: ", round(intercept_val, 4))

  if (intercept_val < 0.8 || intercept_val > 1.5) {
    warning(
      "Unusual LDSC intercept value: ", round(intercept_val, 4),
      "\nTypical values are between 0.8 and 1.5"
    )
  }

  ## ---- Step 4: Apply Z score correction ----
  message("Step 4 of 4: Applying Z score correction")

  n_missing_tstat <- sum(is.na(a$T_STAT))
  if (n_missing_tstat > 0) {
    warning("Found ", n_missing_tstat, " SNPs with missing T_STAT values")
    a <- a[!is.na(a$T_STAT), ]
  }

  if (nrow(a) == 0) {
    stop("No valid SNPs remaining after filtering missing T_STAT values")
  }

  a$z_original <- a$T_STAT
  a$z_adj_int <- a$T_STAT / sqrt(intercept_val)

  a$p_original <- a$P
  a$p_value_int_adj <- 2 * pnorm(-abs(a$z_adj_int))

  n_sig_original <- sum(a$p_original < 0.05, na.rm = TRUE)
  n_sig_adjusted <- sum(a$p_value_int_adj < 0.05, na.rm = TRUE)

  message("  SNPs with p < 0.05 (original): ", n_sig_original)
  message("  SNPs with p < 0.05 (adjusted): ", n_sig_adjusted)

  if (n_sig_original > 0 && n_sig_adjusted < n_sig_original) {
    message(
      "  Correction reduced significant SNPs by ",
      n_sig_original - n_sig_adjusted, " (",
      round(100 * (n_sig_original - n_sig_adjusted) / n_sig_original, 1), "%)"
    )
  }

  ## ---- Save output ----
  output_path <- file.path(tmp_dir, out_file)

  write.table(
    a,
    file = output_path,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  message("Z score correction completed successfully")
  message("  Total SNPs: ", nrow(a))
  message("  Intercept used: ", round(intercept_val, 4))
  message("  Correction factor: ", round(sqrt(intercept_val), 4))
  message("  Results saved to: ", output_path)

  message("  Z score statistics:")
  message("    Original Z mean: ", round(mean(a$z_original, na.rm = TRUE), 4))
  message("    Original Z SD: ", round(sd(a$z_original, na.rm = TRUE), 4))
  message("    Adjusted Z mean: ", round(mean(a$z_adj_int, na.rm = TRUE), 4))
  message("    Adjusted Z SD: ", round(sd(a$z_adj_int, na.rm = TRUE), 4))

  result <- list(
    adjusted_data = a,
    output_file = output_path,
    tmp_dir = tmp_dir,
    intercept_used = intercept_val,
    correction_factor = sqrt(intercept_val),
    n_snps = nrow(a),
    n_sig_original = n_sig_original,
    n_sig_adjusted = n_sig_adjusted
  )

  invisible(result)
}