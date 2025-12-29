#' Perform GWEIS for a binary phenotype and generate a GCIM input file
#'
#' Runs PLINK2 logistic regression with SNP x covariate interaction and exports
#' interaction results (ADDxCOVARk) in GCIM format.
#'
#' @importFrom utils read.table write.table
#' @importFrom stats complete.cases
#' @param plink_path Character. Path to PLINK2 executable.
#' @param tar_mydata Character. Prefix of PLINK binary files (.bed/.bim/.fam).
#' @param tar_pheno_file Character OR list. Phenotype file with columns: FID IID PHENO
#'   (binary PHENO in column 3; typically 0/1 or 1/2). Can pass output from previous step.
#' @param tar_covar_file Character OR list. Covariate file with columns: FID IID COVAR1...
#'   (COVAR1 is column 3, should be PRS). Can pass output from `replace_covariate_with_prs()`.
#' @param int_covar_index Integer. Covariate index for interaction term:
#'   1 -> ADDxCOVAR1 (typically PRS), 2 -> ADDxCOVAR2, ...
#' @param out_file Character. Name of GCIM interaction output file (default "gcim_prs_b_int.txt").
#' @param out_prefix Character. PLINK output prefix inside tmp_dir (default "b_gweis").
#' @param threads Integer. Optional threads for PLINK (default 40).
#'
#' @return A list with file paths and metadata.
#' @export
b_gweis <- function(plink_path,
                    tar_mydata,
                    tar_pheno_file,
                    tar_covar_file,
                    int_covar_index = 1,
                    out_file = "gcim_prs_b_int.txt",
                    out_prefix = "b_gweis",
                    threads = 40) {

  # ---- Helper functions ----
  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Resolve file paths if lists are provided ----
  pheno_source <- NULL
  covar_source <- NULL
  
  if (is.list(tar_pheno_file)) {
    if (!is.null(tar_pheno_file$meta)) {
      pheno_source <- tar_pheno_file$meta$step
    }
    if (!is.null(tar_pheno_file$tar_pheno_file)) {
      tar_pheno_file <- tar_pheno_file$tar_pheno_file
    } else if (!is.null(tar_pheno_file$output_file)) {
      tar_pheno_file <- tar_pheno_file$output_file
    } else {
      .stop2("Input list does not contain 'tar_pheno_file' or 'output_file'.")
    }
  }
  
  if (is.list(tar_covar_file)) {
    if (!is.null(tar_covar_file$meta)) {
      covar_source <- tar_covar_file$meta$step
    }
    if (!is.null(tar_covar_file$output_file)) {
      tar_covar_file <- tar_covar_file$output_file
    } else if (!is.null(tar_covar_file$tar_covar_file)) {
      tar_covar_file <- tar_covar_file$tar_covar_file
    } else {
      .stop2("Input list does not contain 'output_file' or 'tar_covar_file'.")
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

  if (!is.character(tar_pheno_file) || length(tar_pheno_file) != 1 || !nzchar(tar_pheno_file)) {
    .stop2("`tar_pheno_file` must be a non-empty character scalar path.")
  }
  if (!file.exists(tar_pheno_file)) {
    .stop2("Phenotype file not found: ", tar_pheno_file)
  }

  if (!is.character(tar_covar_file) || length(tar_covar_file) != 1 || !nzchar(tar_covar_file)) {
    .stop2("`tar_covar_file` must be a non-empty character scalar path.")
  }
  if (!file.exists(tar_covar_file)) {
    .stop2("Covariate file not found: ", tar_covar_file)
  }

  if (!is.numeric(int_covar_index) || length(int_covar_index) != 1 ||
      is.na(int_covar_index) || int_covar_index < 1) {
    .stop2("`int_covar_index` must be a positive integer")
  }
  int_covar_index <- as.integer(int_covar_index)

  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 1) {
    warning("Invalid `threads` value. Setting to 1.", call. = FALSE)
    threads <- 1
  }
  threads <- as.integer(threads)

  # ---- Check covariate file header and identify covariate name ----
  cov_header <- tryCatch(
    utils::read.table(
      tar_covar_file,
      header = TRUE,
      nrows = 1,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read covariate file header: ", tar_covar_file, "\n", e$message)
  )

  if (!all(c("FID", "IID") %in% names(cov_header))) {
    .stop2(
      "Covariate file must contain 'FID' and 'IID' columns.",
      "\nActual columns: ", paste(names(cov_header), collapse = ", "),
      "\nFile: ", tar_covar_file
    )
  }

  covar_names <- setdiff(names(cov_header), c("FID", "IID"))
  n_covars <- length(covar_names)
  
  if (n_covars < 1) {
    .stop2(
      "Covariate file must contain at least one covariate after FID and IID.",
      "\nFile: ", tar_covar_file
    )
  }
  
  if (int_covar_index > n_covars) {
    .stop2(
      "`int_covar_index` (", int_covar_index, ") exceeds number of covariates (", n_covars, ").",
      "\nAvailable covariates: ", paste(covar_names, collapse = ", ")
    )
  }
  
  covar_name_used <- covar_names[int_covar_index]

  # ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("b_gweis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  int_term <- paste0("ADDxCOVAR", int_covar_index)

  message(.bar())
  message("GCIM-GWEIS-Z: Step 4 - Binary GWEIS (Target Sample)")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Interaction term: ", int_term)
  message("Covariate used: '", covar_name_used, "' (index ", int_covar_index, ")")
  if (!is.null(pheno_source)) message("Phenotype source: ", pheno_source)
  if (!is.null(covar_source)) message("Covariate source: ", covar_source)
  message("Running PLINK2 GWEIS (logistic interaction)...")

  # ---- Run PLINK2 logistic GWEIS ----
  args <- c(
    "--bfile", tar_mydata,
    "--pheno", tar_pheno_file,
    "--pheno-col-nums", "3",
    "--covar", tar_covar_file,
    "--covar-col-nums", "3-",
    "--glm", "interaction", "hide-covar",
    "--allow-no-sex",
    "--covar-variance-standardize",
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

  # ---- Locate GWEIS output ----
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.logistic(\\..*)?$", full.names = TRUE)
  if (length(glm_files) == 0) {
    # Fallback
    glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  }
  if (length(glm_files) == 0) {
    .stop2(
      "No .glm.logistic* output found in: ", tmp_dir,
      "\nCheck PLINK log at: ", log_file
    )
  }
  
  pheno1 <- glm_files[grepl("PHENO1", basename(glm_files))]
  glm_file <- if (length(pheno1) > 0) pheno1[1] else glm_files[1]

  if (!file.exists(glm_file)) {
    .stop2("Selected GWEIS file does not exist: ", glm_file)
  }

  message("GWEIS results file: ", basename(glm_file))

  # ---- Read GWEIS results ----
  gweis_res <- tryCatch(
    utils::read.table(
      glm_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = "",
      check.names = FALSE
    ),
    error = function(e) .stop2("Failed to read GWEIS results from ", glm_file, "\n", e$message)
  )

  if (nrow(gweis_res) == 0) {
    .stop2("GWEIS file has 0 rows: ", glm_file)
  }

  if (!"TEST" %in% names(gweis_res)) {
    .stop2("Column 'TEST' not found in GWEIS output: ", glm_file)
  }

  # ---- Extract interaction effects ----
  gw <- gweis_res[gweis_res$TEST == int_term, , drop = FALSE]
  
  if (nrow(gw) == 0) {
    available_tests <- unique(gweis_res$TEST)
    interaction_tests <- grep("^ADDxCOVAR", available_tests, value = TRUE)
    .stop2(
      "No rows found for interaction term: ", int_term, "\n",
      "Available TEST values: ", paste(available_tests, collapse = ", "), "\n",
      "Available interaction terms: ",
      if (length(interaction_tests) > 0) paste(interaction_tests, collapse = ", ") else "none"
    )
  }

  message("Interaction SNPs extracted: ", nrow(gw))

  # ---- Check required columns ----
  core_needed <- c("ID", "#CHROM", "POS", "A1", "P", "OBS_CT")
  missing_core <- setdiff(core_needed, names(gw))
  if (length(missing_core) > 0) {
    .stop2(
      "Missing required columns in GWEIS results: ",
      paste(missing_core, collapse = ", "),
      "\nColumns available: ", paste(names(gw), collapse = ", ")
    )
  }

  # ---- Identify allele 2 column ----
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(gw)) {
      a2_col <- cand
      break
    }
  }
  if (is.null(a2_col)) {
    .stop2(
      "Could not find a second allele column (tried A2/REF/ALT/ALT1).",
      "\nColumns available: ", paste(names(gw), collapse = ", ")
    )
  }

  message("Using allele 2 column: ", a2_col)

  # ---- Derive BETA: prefer BETA column, else use log(OR) ----
  effect_source <- NULL
  beta <- NULL
  n_removed_effect <- 0
  
  if ("BETA" %in% names(gw)) {
    effect_source <- "BETA"
    beta <- suppressWarnings(as.numeric(gw$BETA))
    bad <- is.na(beta) | is.infinite(beta)
    n_removed_effect <- sum(bad)
    if (any(bad)) {
      warning("Removing ", n_removed_effect, " SNPs with missing/infinite BETA values", call. = FALSE)
      gw <- gw[!bad, , drop = FALSE]
      beta <- beta[!bad]
    }
    
  } else if ("OR" %in% names(gw)) {
    effect_source <- "log(OR)"
    or_num <- suppressWarnings(as.numeric(gw$OR))
    bad <- is.na(or_num) | is.infinite(or_num) | or_num <= 0
    n_removed_or <- sum(bad)
    if (any(bad)) {
      message("Removing ", n_removed_or, " SNPs with missing/invalid OR values")
      gw <- gw[!bad, , drop = FALSE]
      or_num <- or_num[!bad]
    }
    
    # Compute log(OR)
    beta <- log(or_num)
    bad2 <- is.na(beta) | is.infinite(beta)
    n_removed_log <- sum(bad2)
    if (any(bad2)) {
      message("Removing ", n_removed_log, " SNPs with invalid log(OR) values")
      gw <- gw[!bad2, , drop = FALSE]
      beta <- beta[!bad2]
    }
    n_removed_effect <- n_removed_or + n_removed_log
    
  } else {
    .stop2(
      "Neither 'BETA' nor 'OR' column found in logistic GWEIS output.",
      "\nColumns available: ", paste(names(gw), collapse = ", ")
    )
  }

  if (nrow(gw) == 0) {
    .stop2("No valid interaction results after filtering invalid effects.")
  }

  message("Using effect measure: ", effect_source)

  # ---- Build GCIM data frame ----
  gcim_data <- data.frame(
    order = seq_len(nrow(gw)),
    snpid = gw$ID,
    chr   = gw$`#CHROM`,
    bp    = gw$POS,
    a1    = gw$A1,
    a2    = gw[[a2_col]],
    beta  = beta,
    pval  = gw$P,
    N     = gw$OBS_CT,
    stringsAsFactors = FALSE
  )

  # ---- Remove incomplete cases ----
  keep <- stats::complete.cases(gcim_data)
  n_removed_missing <- sum(!keep)
  
  if (!all(keep)) {
    warning("Removing ", n_removed_missing, " SNPs with missing values in key columns", call. = FALSE)
    gcim_data <- gcim_data[keep, , drop = FALSE]
  }
  
  if (nrow(gcim_data) == 0) {
    .stop2("No valid interaction results after filtering missing values.")
  }

  # ---- Write GCIM file ----
  gcim_file <- file.path(tmp_dir, out_file)
  writeLines("order snpid chr bp a1 a2 beta pval N", con = gcim_file)
  utils::write.table(
    gcim_data,
    file = gcim_file,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = " ",
    append = TRUE
  )

  if (!file.exists(gcim_file)) {
    .stop2("Failed to create GCIM file: ", gcim_file)
  }

  n_removed_total <- n_removed_effect + n_removed_missing

  message(.bar())
  message("Binary GWEIS completed successfully!")
  message("  GLM file: ", glm_file)
  message("  Interaction term: ", int_term)
  message("  Covariate: '", covar_name_used, "'")
  message("  Effect measure: ", effect_source)
  message("  Interaction SNPs: ", nrow(gcim_data))
  if (n_removed_total > 0) {
    message("  SNPs removed (invalid effects): ", n_removed_effect)
    message("  SNPs removed (missing data): ", n_removed_missing)
    message("  Total SNPs removed: ", n_removed_total)
  }
  message("  GCIM file: ", gcim_file)
  message(.bar())

  # ---- Return results (pipeline-friendly structure) ----
  res <- list(
    # Backwards compatible fields
    glm_file = glm_file,
    gcim_file = gcim_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_interactions = nrow(gcim_data),
    interaction_term = int_term,
    covariate_name = covar_name_used,

    # Pipeline-friendly additions
    files = list(
      glm = glm_file,
      gcim = gcim_file,
      plink_log = log_file,
      stdout = stdout_file,
      stderr = stderr_file
    ),
    meta = list(
      step = "GWEIS",
      trait_type = "binary",
      model = "logistic_interaction",
      int_covar_index = int_covar_index,
      int_term = int_term,
      covariate_name = covar_name_used,
      n_covariates = n_covars,
      threads = threads,
      a2_col = a2_col,
      effect_source = effect_source,
      n_snps_total = nrow(gweis_res),
      n_snps_interaction = nrow(gcim_data),
      n_snps_removed_effect = n_removed_effect,
      n_snps_removed_missing = n_removed_missing,
      n_snps_removed_total = n_removed_total,
      pheno_source = pheno_source,
      covar_source = covar_source,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}