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
#' @param int_term Character OR NULL. Interaction term to extract (e.g., "ADDxCOVAR1").
#'   If NULL and a list is provided, attempts to extract from metadata.
#' @param out_file Character. Output filename for adjusted full table (default "gcim_z_adjusted.txt").
#' @param out_ldsc_name Character. Output filename for minimal adjusted LDSC-ready table
#'   (default "ldsc_ready_z_adjusted.tsv").
#'
#' @return A list with adjusted results and output paths.
#' @export
gcim_z_adjust <- function(glm_file,
                          intercept_file,
                          int_term = NULL,
                          out_file = "gcim_z_adjusted.txt",
                          out_ldsc_name = "ldsc_ready_z_adjusted.tsv") {

  # ---- Helper functions ----
  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Resolve inputs if lists are provided ----
  gweis_source <- NULL
  ldsc_source <- NULL
  
  # Extract GLM file path
  if (is.list(glm_file)) {
    if (!is.null(glm_file$meta)) {
      gweis_source <- glm_file$meta$trait_type
      # Auto-extract interaction term if not provided
      if (is.null(int_term) && !is.null(glm_file$meta$int_term)) {
        int_term <- glm_file$meta$int_term
        message("Auto-detected interaction term: ", int_term)
      }
    }
    
    if (!is.null(glm_file$glm_file)) {
      glm_file <- glm_file$glm_file
    } else if (!is.null(glm_file$files$glm)) {
      glm_file <- glm_file$files$glm
    } else {
      .stop2("Input list does not contain 'glm_file' or 'files$glm' element.")
    }
  }

  # Extract intercept value
  intercept_val_from_list <- NULL
  intercept_se <- NULL
  
  if (is.list(intercept_file)) {
    if (!is.null(intercept_file$meta)) {
      ldsc_source <- intercept_file$meta$step
    }
    
    if (!is.null(intercept_file$intercept) && is.finite(intercept_file$intercept)) {
      intercept_val_from_list <- as.numeric(intercept_file$intercept)
      if (!is.null(intercept_file$intercept_se) && is.finite(intercept_file$intercept_se)) {
        intercept_se <- as.numeric(intercept_file$intercept_se)
      }
    }
    
    if (!is.null(intercept_file$intercept_file)) {
      intercept_file <- intercept_file$intercept_file
    } else if (!is.null(intercept_file$files$intercept_txt)) {
      intercept_file <- intercept_file$files$intercept_txt
    }
  }

  # ---- Input validation ----
  if (!is.character(glm_file) || length(glm_file) != 1 || !nzchar(glm_file)) {
    .stop2("`glm_file` must be a non-empty character scalar path.")
  }
  if (!file.exists(glm_file)) {
    .stop2("GWEIS GLM file not found: ", glm_file)
  }

  if (is.null(int_term) || !is.character(int_term) || 
      length(int_term) != 1 || !nzchar(int_term)) {
    .stop2(
      "`int_term` must be a non-empty string (e.g., 'ADDxCOVAR1').",
      "\nIf passing a list from q_gweis/b_gweis, the term should be auto-extracted."
    )
  }

  # Parse intercept value
  if (!is.null(intercept_val_from_list)) {
    intercept_val <- intercept_val_from_list
  } else {
    if (!is.character(intercept_file) || length(intercept_file) != 1 || !nzchar(intercept_file)) {
      .stop2("`intercept_file` must be a non-empty character scalar path.")
    }
    if (!file.exists(intercept_file)) {
      .stop2("LDSC intercept file not found: ", intercept_file)
    }
    
    intercept_lines <- tryCatch(
      readLines(intercept_file, warn = FALSE),
      error = function(e) .stop2("Failed to read intercept file: ", intercept_file, "\n", e$message)
    )
    
    intercept_val <- suppressWarnings(as.numeric(trimws(intercept_lines[1])))

    # Fallback if file format is "Intercept: 1.02"
    if (is.na(intercept_val)) {
      idx <- grep("Intercept:", intercept_lines, ignore.case = TRUE)
      if (length(idx) > 0) {
        line <- intercept_lines[idx[1]]
        intercept_val <- suppressWarnings(
          as.numeric(sub("^.*Intercept:\\s*([0-9.eE+-]+).*$", "\\1", line))
        )
      }
    }
    
    # Try to extract SE if present
    se_idx <- grep("^SE:", intercept_lines, ignore.case = TRUE)
    if (length(se_idx) > 0) {
      intercept_se <- suppressWarnings(
        as.numeric(sub("^.*SE:\\s*([0-9.eE+-]+).*$", "\\1", intercept_lines[se_idx[1]]))
      )
    }
  }

  if (!is.finite(intercept_val) || intercept_val <= 0) {
    .stop2(
      "Invalid LDSC intercept value: ", intercept_val,
      "\nIntercept must be a positive finite number."
    )
  }

  # ---- Create temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("gcim_z_adjust_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  message(.bar())
  message("GCIM-GWEIS-Z: Step 7 - Z-Score Adjustment")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Interaction term: ", int_term)
  if (!is.null(gweis_source)) {
    message("GWEIS source: ", gweis_source, " trait")
  }
  if (!is.null(ldsc_source)) {
    message("LDSC source: ", ldsc_source)
  }
  message("LDSC intercept: ", sprintf("%.6f", intercept_val))
  if (!is.null(intercept_se) && is.finite(intercept_se)) {
    message("Intercept SE: ", sprintf("%.6f", intercept_se))
  }

  # ---- Read GWEIS results ----
  gweis_res <- tryCatch(
    utils::read.table(
      glm_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = "",
      check.names = FALSE
    ),
    error = function(e) .stop2("Failed to read GWEIS file: ", glm_file, "\n", e$message)
  )

  if (nrow(gweis_res) == 0) {
    .stop2("GWEIS file has 0 rows: ", glm_file)
  }

  if (!"TEST" %in% names(gweis_res)) {
    .stop2("Column 'TEST' not found in GWEIS results: ", glm_file)
  }

  # ---- Extract interaction term ----
  a <- gweis_res[gweis_res$TEST == int_term, , drop = FALSE]
  
  if (nrow(a) == 0) {
    available_tests <- unique(gweis_res$TEST)
    interaction_tests <- grep("^ADDxCOVAR", available_tests, value = TRUE)
    .stop2(
      "No interaction term found: ", int_term,
      "\nAvailable TEST values: ", paste(utils::head(available_tests, 20), collapse = ", "),
      "\nAvailable interaction terms: ",
      if (length(interaction_tests) > 0) paste(interaction_tests, collapse = ", ") else "none"
    )
  }

  message("Interaction SNPs: ", nrow(a))

  # ---- Determine Z-score source ----
  z_source <- NULL
  
  if ("Z_STAT" %in% names(a)) {
    z_source <- "Z_STAT"
    z <- suppressWarnings(as.numeric(a$Z_STAT))
  } else if ("T_STAT" %in% names(a)) {
    z_source <- "T_STAT"
    z <- suppressWarnings(as.numeric(a$T_STAT))
  } else if (all(c("BETA", "SE") %in% names(a))) {
    z_source <- "BETA/SE"
    beta <- suppressWarnings(as.numeric(a$BETA))
    se <- suppressWarnings(as.numeric(a$SE))
    z <- beta / se
  } else {
    .stop2(
      "Cannot compute Z: need Z_STAT or T_STAT or (BETA and SE).",
      "\nColumns available: ", paste(names(a), collapse = ", ")
    )
  }

  message("Z-score source: ", z_source)

  # ---- Extract original p-values ----
  p_orig <- if ("P" %in% names(a)) {
    suppressWarnings(as.numeric(a$P))
  } else {
    2 * stats::pnorm(-abs(z))
  }

  # ---- Filter invalid statistics ----
  keep <- is.finite(z) & is.finite(p_orig)
  n_removed_invalid <- sum(!keep)
  
  if (!all(keep)) {
    warning(
      "Dropping ", n_removed_invalid, " SNPs with missing/invalid Z or P values.",
      call. = FALSE
    )
    a <- a[keep, , drop = FALSE]
    z <- z[keep]
    p_orig <- p_orig[keep]
  }
  
  if (nrow(a) == 0) {
    .stop2("No valid SNPs remaining after filtering invalid statistics.")
  }

  # ---- Remove duplicate SNP IDs ----
  n_duplicates <- 0
  if ("ID" %in% names(a) && anyDuplicated(a$ID)) {
    n_duplicates <- sum(duplicated(a$ID))
    warning(
      "Found ", n_duplicates, " duplicated SNP IDs; keeping first occurrence.",
      call. = FALSE
    )
    keep2 <- !duplicated(a$ID)
    a <- a[keep2, , drop = FALSE]
    z <- z[keep2]
    p_orig <- p_orig[keep2]
  }

  # ---- Apply LDSC intercept correction ----
  corr_factor <- sqrt(intercept_val)
  z_adj <- z / corr_factor
  p_adj <- 2 * stats::pnorm(-abs(z_adj))

  # Add new columns
  a$z_original <- z
  a$z_adj_int <- z_adj
  a$p_original <- p_orig
  a$p_value_int_adj <- p_adj
  a$ldsc_intercept_used <- intercept_val

  # ---- Compute significance statistics ----
  n_sig_original <- sum(a$p_original < 0.05, na.rm = TRUE)
  n_sig_adjusted <- sum(a$p_value_int_adj < 0.05, na.rm = TRUE)
  n_sig_bonf_orig <- sum(a$p_original < 0.05 / nrow(a), na.rm = TRUE)
  n_sig_bonf_adj <- sum(a$p_value_int_adj < 0.05 / nrow(a), na.rm = TRUE)

  # ---- Save full adjusted output ----
  output_path <- file.path(tmp_dir, out_file)
  utils::write.table(
    a,
    file = output_path,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  if (!file.exists(output_path)) {
    .stop2("Failed to create adjusted output file: ", output_path)
  }

  # ---- Create minimal LDSC-ready adjusted table ----
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(a)) {
      a2_col <- cand
      break
    }
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
      utils::write.table(
        ldsc_adj,
        file = ldsc_adj_path,
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
      )
    } else {
      ldsc_adj_path <- NA_character_
    }
  } else {
    ldsc_adj_path <- NA_character_
  }

  message(.bar())
  message("Z-score adjustment completed successfully!")
  message("  Correction factor (sqrt(intercept)): ", sprintf("%.6f", corr_factor))
  message("  Valid SNPs: ", nrow(a))
  if (n_removed_invalid > 0) {
    message("  SNPs removed (invalid): ", n_removed_invalid)
  }
  if (n_duplicates > 0) {
    message("  SNPs removed (duplicates): ", n_duplicates)
  }
  message("")
  message("Significance counts (p < 0.05):")
  message("  Original:  ", n_sig_original, " / ", nrow(a), 
          " (", sprintf("%.2f%%", 100 * n_sig_original / nrow(a)), ")")
  message("  Adjusted:  ", n_sig_adjusted, " / ", nrow(a), 
          " (", sprintf("%.2f%%", 100 * n_sig_adjusted / nrow(a)), ")")
  message("")
  message("Bonferroni-corrected (p < ", sprintf("%.2e", 0.05 / nrow(a)), "):")
  message("  Original:  ", n_sig_bonf_orig)
  message("  Adjusted:  ", n_sig_bonf_adj)
  message("")
  message("  Full results: ", output_path)
  if (!is.na(ldsc_adj_path)) {
    message("  LDSC-ready (adjusted): ", ldsc_adj_path)
  }
  message(.bar())

  # ---- Return results (pipeline-friendly structure) ----
  res <- list(
    # Backwards compatible fields
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

    # Pipeline-friendly additions
    files = list(
      adjusted_full = output_path,
      ldsc_ready_adjusted = ldsc_adj_path
    ),
    meta = list(
      step = "Z_ADJUST",
      int_term = int_term,
      gweis_source = gweis_source,
      ldsc_source = ldsc_source,
      intercept = intercept_val,
      intercept_se = if (!is.null(intercept_se) && is.finite(intercept_se)) intercept_se else NA_real_,
      corr_factor = corr_factor,
      z_source = z_source,
      n_snps_input = nrow(a) + n_removed_invalid + n_duplicates,
      n_snps_valid = nrow(a),
      n_snps_removed_invalid = n_removed_invalid,
      n_snps_removed_duplicates = n_duplicates,
      n_sig_original = n_sig_original,
      n_sig_adjusted = n_sig_adjusted,
      n_sig_bonf_original = n_sig_bonf_orig,
      n_sig_bonf_adjusted = n_sig_bonf_adj,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}