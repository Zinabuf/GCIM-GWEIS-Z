#' Perform GWEIS for a binary phenotype and generate a GCIM input file
#'
#' Runs PLINK2 logistic regression with SNP x covariate interaction and exports
#' interaction results (ADDxCOVARk) in GCIM format.
#'
#' @importFrom utils read.table write.table
#' @importFrom stats complete.cases
#' @param plink_path Character. Path to PLINK2 executable.
#' @param tar_mydata Character. Prefix of PLINK binary files (.bed/.bim/.fam).
#' @param tar_pheno_file Character. Phenotype file with columns: FID IID PHENO
#'   (binary PHENO in column 3; typically 0/1 or 1/2).
#' @param tar_covar_file Character. Covariate file with columns: FID IID COVAR1...
#'   (COVAR1 is column 3, used as COVAR1 in PLINK naming).
#' @param int_covar_index Integer. Covariate index for interaction term:
#'   1 -> ADDxCOVAR1, 2 -> ADDxCOVAR2, ...
#' @param out_file Character. Name of GCIM interaction output file.
#' @param out_prefix Character. PLINK output prefix inside tmp_dir (default "b_gweis").
#' @param threads Integer. Optional threads for PLINK (default 1).
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
                    threads = 1) {

  ## ---- Input validation ----
  if (!file.exists(plink_path)) stop("PLINK executable not found at: ", plink_path)

  for (ext in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(tar_mydata, ext))) {
      stop("Missing PLINK file: ", paste0(tar_mydata, ext))
    }
  }

  if (!file.exists(tar_pheno_file)) stop("Phenotype file not found: ", tar_pheno_file)
  if (!file.exists(tar_covar_file)) stop("Covariate file not found: ", tar_covar_file)

  if (!is.numeric(int_covar_index) || length(int_covar_index) != 1 ||
      is.na(int_covar_index) || int_covar_index < 1) {
    stop("int_covar_index must be a positive integer")
  }

  ## ---- Check covariate header / count ----
  cov_header <- tryCatch(
    utils::read.table(tar_covar_file,
               header = TRUE,
               nrows = 1,
               stringsAsFactors = FALSE,
               check.names = FALSE),
    error = function(e) stop("Failed to read covariate file header: ", tar_covar_file, "\n", e$message)
  )

  if (!all(c("FID", "IID") %in% names(cov_header))) {
    stop("Covariate file must contain 'FID' and 'IID' columns")
  }

  n_covars <- ncol(cov_header) - 2
  if (n_covars < 1) stop("Covariate file must contain at least one covariate after FID and IID")
  if (int_covar_index > n_covars) {
    stop("int_covar_index (", int_covar_index, ") exceeds number of covariates (", n_covars, ")")
  }

  ## ---- Temp dir ----
  tmp_dir <- file.path(tempdir(), paste0("b_gweis_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", tmp_dir)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")

  ## ---- Run PLINK2 logistic GWEIS ----
  int_term <- paste0("ADDxCOVAR", int_covar_index)
  message("Step 1 of 3: Running PLINK logistic GWEIS (interaction term: ", int_term, ")")

  args <- c(
    "--bfile", tar_mydata,
    "--pheno", tar_pheno_file,
    "--pheno-col-nums", "3",
    "--covar", tar_covar_file,
    "--covar-col-nums", "3-n",
    "--glm", "interaction", "hide-covar",
    "--allow-no-sex",
    "--covar-variance-standardize",
    # Phenotype variance standardization is not applicable to logistic; harmless if ignored in some versions,
    # but we omit it to avoid confusion.
    "--threads", as.character(threads),
    "--out", out_pref
  )

  exit_code <- suppressWarnings(system2(plink_path, args = args))
  if (!is.null(exit_code) && exit_code != 0) {
    warning("PLINK exited with non-zero status (", exit_code, "). Check: ", log_file)
  }

  ## ---- Locate GWEIS output robustly ----
  message("Step 2 of 3: Locating logistic GWEIS output")
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.logistic(\\..*)?$", full.names = TRUE)
  if (length(glm_files) == 0) {
    stop("No .glm.logistic* output found in: ", tmp_dir, "\nCheck PLINK log at: ", log_file)
  }
  glm_file <- glm_files[grepl("PHENO1", basename(glm_files))]
  if (length(glm_file) == 0) glm_file <- glm_files[1]
  glm_file <- glm_file[1]

  ## ---- Read results ----
  gweis_res <- tryCatch(
    utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE),
    error = function(e) stop("Failed to read GWEIS results from ", glm_file, "\n", e$message)
  )

  if (!"TEST" %in% names(gweis_res)) stop("Column 'TEST' not found in GWEIS output: ", glm_file)

  ## ---- Extract interaction effects ----
  message("Step 3 of 3: Extracting interaction effects (", int_term, ")")

  tests <- unique(gweis_res$TEST)
  gw <- gweis_res[gweis_res$TEST == int_term, , drop = FALSE]

  if (nrow(gw) == 0) {
    stop("No rows found for interaction term: ", int_term, "\nAvailable interaction terms: ",
         paste(grep("^ADDxCOVAR", tests, value = TRUE), collapse = ", "))
  }

  ## ---- Column mapping ----
  # We need: ID, #CHROM, POS, A1, P, OBS_CT, and an effect size.
  core_needed <- c("ID", "#CHROM", "POS", "A1", "P", "OBS_CT")
  missing_core <- setdiff(core_needed, names(gw))
  if (length(missing_core) > 0) {
    stop("Missing required columns in GWEIS results: ",
         paste(missing_core, collapse = ", "),
         "\nColumns available: ", paste(names(gw), collapse = ", "))
  }

  # Determine allele2 column
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(gw)) { a2_col <- cand; break }
  }
  if (is.null(a2_col)) {
    stop("Could not find a second allele column (tried A2/REF/ALT/ALT1). Columns: ",
         paste(names(gw), collapse = ", "))
  }

  # Determine effect: prefer BETA; else use log(OR)
  if ("BETA" %in% names(gw)) {
    beta <- suppressWarnings(as.numeric(gw$BETA))
    bad <- is.na(beta) | is.infinite(beta)
    if (any(bad)) {
      warning("Removing ", sum(bad), " SNPs with invalid BETA")
      gw <- gw[!bad, , drop = FALSE]
      beta <- beta[!bad]
    }
  } else if ("OR" %in% names(gw)) {
    or_num <- suppressWarnings(as.numeric(gw$OR))
    bad <- is.na(or_num) | or_num <= 0
    if (any(bad)) {
      warning("Removing ", sum(bad), " SNPs with missing/non-positive OR")
      gw <- gw[!bad, , drop = FALSE]
      or_num <- or_num[!bad]
    }
    beta <- log(or_num)
    bad2 <- is.na(beta) | is.infinite(beta)
    if (any(bad2)) {
      warning("Removing ", sum(bad2), " SNPs with invalid log(OR)")
      gw <- gw[!bad2, , drop = FALSE]
      beta <- beta[!bad2]
    }
  } else {
    stop("Neither 'BETA' nor 'OR' found in logistic GWEIS output. Columns: ",
         paste(names(gw), collapse = ", "))
  }

  if (nrow(gw) == 0) stop("No valid interaction results after filtering invalid effects.")

  ## ---- Build GCIM data ----
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

  # Drop incomplete
  keep <- stats::complete.cases(gcim_data)
  if (!all(keep)) {
    warning("Removing ", sum(!keep), " SNPs with missing values")
    gcim_data <- gcim_data[keep, , drop = FALSE]
  }
  if (nrow(gcim_data) == 0) stop("No valid interaction results after filtering missing values.")

  ## ---- Write GCIM file ----
  gcim_file <- file.path(tmp_dir, out_file)
  writeLines("order snpid chr bp a1 a2 beta pval N", con = gcim_file)
  utils::write.table(gcim_data, file = gcim_file, quote = FALSE, col.names = FALSE,
              row.names = FALSE, sep = " ", append = TRUE)

  message("Binary GWEIS completed")
  message("  GLM file: ", glm_file)
  message("  Interaction term: ", int_term)
  message("  Interaction SNPs: ", nrow(gcim_data))
  message("  GCIM file: ", gcim_file)

  invisible(list(
    glm_file = glm_file,
    gcim_file = gcim_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_interactions = nrow(gcim_data),
    interaction_term = int_term
  ))
}