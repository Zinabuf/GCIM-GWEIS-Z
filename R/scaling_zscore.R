#' LDSC intercept based Z score correction for GWEIS interaction results
#'
#' Applies LDSC intercept correction to interaction test statistics:
#'   Z_adj = Z_original / sqrt(intercept)
#' and recomputes two-sided p-values.
#'
#' @importFrom utils read.table write.table
#' @param glm_file Character OR list. PLINK GWEIS output file (.glm.linear or .glm.logistic*)
#'   OR list from `q_gweis()`/`b_gweis()` containing `$glm_file` or `$files$glm`.
#' @param intercept_file Character OR list. File containing LDSC intercept value from `ldsc_h2_gcim()`,
#'   OR list from `ldsc_h2_gcim()` containing `$intercept` / `$intercept_file` / `$files$intercept_txt`.
#' @param int_term Character. Interaction term to extract (e.g., "ADDxCOVAR1").
#' @param out_file Character. Output filename for adjusted full table (default "gcim_z_adjusted.txt").
#' @param out_ldsc_name Character. Output filename for minimal adjusted LDSC-ready table
#'   (default "ldsc_ready_z_adjusted.tsv").
#'
#' @return A list with adjusted results and output paths.
#' @export
gcim_z_adjust <- function(glm_file,
                          intercept_file,
                          int_term = "ADDxCOVAR1",
                          out_file = "gcim_z_adjusted.txt",
                          out_ldsc_name = "ldsc_ready_z_adjusted.tsv") {

  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Allow passing prior-step lists ----
  if (is.list(glm_file) && !is.null(glm_file$glm_file)) glm_file <- glm_file$glm_file
  if (is.list(glm_file) && !is.null(glm_file$files$glm)) glm_file <- glm_file$files$glm

  intercept_val_from_list <- NULL
  if (is.list(intercept_file)) {
    if (!is.null(intercept_file$intercept) && is.finite(intercept_file$intercept)) {
      intercept_val_from_list <- as.numeric(intercept_file$intercept)
    }
    if (!is.null(intercept_file$intercept_file)) intercept_file <- intercept_file$intercept_file
    if (!is.null(intercept_file$files$intercept_txt)) intercept_file <- intercept_file$files$intercept_txt
  }

  # ---- Input validation ----
  if (!is.character(glm_file) || length(glm_file) != 1 || !file.exists(glm_file)) {
    .stop2("GWEIS GLM file not found: ", glm_file)
  }

  if (!is.null(intercept_val_from_list)) {
    intercept_val <- intercept_val_from_list
  } else {
    if (!is.character(intercept_file) || length(intercept_file) != 1 || !file.exists(intercept_file)) {
      .stop2("LDSC intercept file not found: ", intercept_file)
    }
    intercept_lines <- tryCatch(
      readLines(intercept_file, warn = FALSE),
      error = function(e) .stop2("Failed to read intercept file: ", intercept_file, "\n", e$message)
    )
    intercept_val <- suppressWarnings(as.numeric(trimws(intercept_lines[1])))

    # Fallback if file is like "Intercept: 1.02"
    if (is.na(intercept_val)) {
      idx <- grep("Intercept:", intercept_lines, ignore.case = TRUE)
      if (length(idx) > 0) {
        line <- intercept_lines[idx[1]]
        intercept_val <- suppressWarnings(as.numeric(sub("^.*Intercept:\\s*([0-9.eE+-]+).*$", "\\1", line)))
      }
    }
  }

  if (!is.finite(intercept_val) || intercept_val <= 0) {
    .stop2("Invalid LDSC intercept value: ", intercept_val)
  }

  # ---- Temp dir ----
  tmp_dir <- file.path(tempdir(), paste0("gcim_z_adjust_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Read GLM ----
  gweis_res <- tryCatch(
    utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE,
                      comment.char = "", check.names = FALSE),
    error = function(e) .stop2("Failed to read GWEIS file: ", glm_file, "\n", e$message)
  )
  if (!"TEST" %in% names(gweis_res)) .stop2("Column 'TEST' not found in GWEIS results: ", glm_file)

  # ---- Extract interaction term ----
  a <- gweis_res[gweis_res$TEST == int_term, , drop = FALSE]
  if (nrow(a) == 0) {
    available_tests <- unique(gweis_res$TEST)
    .stop2("No interaction term found: ", int_term, "\nAvailable interaction terms: ",
           paste(grep("^ADDxCOVAR", available_tests, value = TRUE), collapse = ", "))
  }

  # ---- Z source ----
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
    .stop2("Cannot compute Z: need Z_STAT or T_STAT or (BETA and SE). Columns: ",
           paste(names(a), collapse = ", "))
  }

  # Original p-values
  p_orig <- if ("P" %in% names(a)) suppressWarnings(as.numeric(a$P)) else 2 * stats::pnorm(-abs(z))

  # Drop invalid stats
  keep <- is.finite(z) & is.finite(p_orig)
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " SNPs with missing/invalid Z or P.", call. = FALSE)
    a <- a[keep, , drop = FALSE]
    z <- z[keep]
    p_orig <- p_orig[keep]
  }
  if (nrow(a) == 0) .stop2("No valid SNPs remaining after filtering invalid statistics.")

  # Drop duplicate SNP IDs (safer downstream)
  if ("ID" %in% names(a) && anyDuplicated(a$ID)) {
    warning("Duplicated SNP IDs detected; keeping first occurrence.", call. = FALSE)
    keep2 <- !duplicated(a$ID)
    a <- a[keep2, , drop = FALSE]
    z <- z[keep2]
    p_orig <- p_orig[keep2]
  }

  # ---- Apply correction ----
  corr_factor <- sqrt(intercept_val)
  z_adj <- z / corr_factor
  p_adj <- 2 * stats::pnorm(-abs(z_adj))

  a$z_original <- z
  a$z_adj_int <- z_adj
  a$p_original <- p_orig
  a$p_value_int_adj <- p_adj
  a$ldsc_intercept_used <- intercept_val

  n_sig_original <- sum(a$p_original < 0.05, na.rm = TRUE)
  n_sig_adjusted <- sum(a$p_value_int_adj < 0.05, na.rm = TRUE)

  # ---- Save full output ----
  output_path <- file.path(tmp_dir, out_file)
  utils::write.table(a, file = output_path, row.names = FALSE, quote = FALSE, sep = "\t")

  # ---- Also write minimal LDSC-ready adjusted table (optional but useful) ----
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(a)) { a2_col <- cand; break }
  }
  ldsc_adj_path <- file.path(tmp_dir, out_ldsc_name)

  if (!is.null(a2_col) && all(c("ID", "A1", "OBS_CT") %in% names(a))) {
    ldsc_adj <- data.frame(
      SNP = as.character(a$ID),
      A1 = as.character(a$A1),
      A2 = as.character(a[[a2_col]]),
      N  = suppressWarnings(as.numeric(a$OBS_CT)),
      Z  = as.numeric(a$z_adj_int),
      stringsAsFactors = FALSE
    )
    keep_ld <- stats::complete.cases(ldsc_adj) & is.finite(ldsc_adj$Z) & ldsc_adj$N > 0
    ldsc_adj <- ldsc_adj[keep_ld, , drop = FALSE]
    if (nrow(ldsc_adj) > 0) {
      utils::write.table(ldsc_adj, file = ldsc_adj_path, row.names = FALSE,
                         col.names = TRUE, sep = "\t", quote = FALSE)
    } else {
      ldsc_adj_path <- NA_character_
    }
  } else {
    ldsc_adj_path <- NA_character_
  }

  invisible(list(
    adjusted_data = a,
    output_file = output_path,
    ldsc_adjusted_sumstats = ldsc_adj_path,
    tmp_dir = tmp_dir,
    intercept_used = intercept_val,
    correction_factor = corr_factor,
    z_source = z_source,
    n_snps = nrow(a),
    n_sig_original = n_sig_original,
    n_sig_adjusted = n_sig_adjusted,
    files = list(
      adjusted_full = output_path,
      ldsc_ready_adjusted = ldsc_adj_path
    ),
    meta = list(
      step = "Z_ADJUST",
      int_term = int_term,
      intercept = intercept_val,
      corr_factor = corr_factor,
      z_source = z_source
    )
  ))
}