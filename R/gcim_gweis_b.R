#' Perform GWEIS for a binary phenotype and generate a GCIM input file
#'
#' This function runs a PLINK logistic GWEIS (SNP by covariate interaction)
#' and outputs a GCIM-formatted interaction file containing SNP position,
#' alleles, interaction effect size on the log odds scale, p value, and
#' sample size. All outputs are saved to a temporary directory.
#'
#' @param plink_path Character. Path to the PLINK (v2.0) executable.
#' @param dis_snp Character. Prefix of PLINK binary files for discovery data.
#' @param bp_dis_phen Character. File path to phenotype file (FID IID PHENO).
#'   Binary phenotype should be coded as 0/1 or 1/2.
#' @param bp_dis_cov Character. File path to covariate file with PRS
#'   (typically output from \code{replace_covariate_with_prs()}).
#'   Format: FID IID COVAR1 COVAR2 ... COVARn
#' @param int_covar_index Integer. Covariate column index for interaction
#'   (for example, 1 for ADDxCOVAR1). Default is 1.
#' @param out_file Character. Name of GCIM interaction output file.
#'   Default is "gcim_prs_b_int.txt".
#'
#' @return A list containing file paths to GWEIS outputs stored in a temporary
#' directory:
#' \describe{
#'   \item{glm_file}{Path to PLINK .glm.logistic.hybrid file}
#'   \item{gcim_file}{Path to GCIM-formatted interaction file}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_interactions}{Number of SNPs with interaction results}
#'   \item{interaction_term}{Name of the interaction term tested}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Validates input files (PLINK binaries, phenotype, covariates)
#'   \item Creates a temporary directory for all outputs
#'   \item Runs PLINK logistic GWEIS with an interaction term
#'   \item Extracts interaction effects (ADDxCOVARn)
#'   \item Converts odds ratios to log odds for BETA values
#'   \item Formats results for GCIM-GWEIS-Z analysis
#'   \item Writes a GCIM-formatted file with a header line
#' }
#'
#' The GCIM output file format is:
#' order snpid chr bp a1 a2 beta pval N
#'
#' @export
b_gweis <- function(plink_path,
                    dis_snp,
                    bp_dis_phen,
                    bp_dis_cov,
                    int_covar_index = 1,
                    out_file = "gcim_prs_b_int.txt") {

  ## ---- Input validation ----
  if (!file.exists(plink_path)) {
    stop("PLINK executable not found at: ", plink_path)
  }

  if (!file.exists(paste0(dis_snp, ".bed")) ||
      !file.exists(paste0(dis_snp, ".bim")) ||
      !file.exists(paste0(dis_snp, ".fam"))) {
    stop("PLINK binary files (.bed/.bim/.fam) not found with prefix: ", dis_snp)
  }

  if (!file.exists(bp_dis_phen)) {
    stop("Phenotype file not found: ", bp_dis_phen)
  }

  if (!file.exists(bp_dis_cov)) {
    stop("Covariate file not found: ", bp_dis_cov)
  }

  ## ---- Validate covariate index ----
  if (!is.numeric(int_covar_index) || length(int_covar_index) != 1 ||
      is.na(int_covar_index) || int_covar_index < 1) {
    stop("int_covar_index must be a positive integer")
  }

  ## ---- Check covariate file has enough columns ----
  cov_header <- tryCatch(
    {
      read.table(bp_dis_cov,
                 header = TRUE,
                 nrows = 1,
                 stringsAsFactors = FALSE,
                 check.names = FALSE)
    },
    error = function(e) {
      stop("Failed to read covariate file header: ", bp_dis_cov, "\n", e$message)
    }
  )

  if (!all(c("FID", "IID") %in% names(cov_header))) {
    stop("Covariate file must contain 'FID' and 'IID' columns")
  }

  n_covars <- ncol(cov_header) - 2
  if (n_covars < 1) {
    stop("Covariate file must contain at least one covariate column after FID and IID")
  }

  if (int_covar_index > n_covars) {
    stop(
      "int_covar_index (", int_covar_index, ") exceeds number of covariates (",
      n_covars, ") in file: ", bp_dis_cov
    )
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("b_gweis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }
  message("Output directory: ", tmp_dir)

  out_prefix <- file.path(tmp_dir, "b_gweis")

  ## ---- Build --parameters list (explicit indices) ----
  param_max <- 2 + n_covars
  param_list <- paste(seq_len(param_max), collapse = ",")

  ## ---- Run PLINK logistic GWEIS ----
  message("Step 1 of 3: Running PLINK logistic GWEIS with interaction term")

  int_term <- paste0("ADDxCOVAR", int_covar_index)
  message("  Interaction term: ", int_term)

  cmd <- paste(
    shQuote(plink_path),
    "--bfile", shQuote(dis_snp),
    "--pheno", shQuote(bp_dis_phen),
    "--glm interaction",
    "--covar", shQuote(bp_dis_cov),
    "--parameters", shQuote(param_list),
    "--allow-no-sex",
    "--covar-variance-standardize",
    "--out", shQuote(out_prefix)
  )

  exit_code <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)

  if (exit_code != 0) {
    warning(
      "PLINK exited with non-zero status. Check log file: ",
      paste0(out_prefix, ".log")
    )
  }

  ## ---- Locate GWEIS results ----
  message("Step 2 of 3: Reading logistic GWEIS results")

  glm_file <- paste0(out_prefix, ".PHENO1.glm.logistic.hybrid")
  log_file <- paste0(out_prefix, ".log")

  if (!file.exists(glm_file)) {
    stop(
      "Binary GWEIS output file not found: ", glm_file,
      "\nCheck PLINK log at: ", log_file
    )
  }

  ## ---- Read GWEIS results ----
  gweis_res <- tryCatch(
    {
      read.table(
        glm_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        comment.char = ""
      )
    },
    error = function(e) {
      stop("Failed to read GWEIS results from ", glm_file, "\n", e$message)
    }
  )

  if (!"TEST" %in% colnames(gweis_res)) {
    stop("Column 'TEST' not found in GWEIS results")
  }

  ## ---- Extract interaction term ----
  message("Step 3 of 3: Extracting and converting interaction effects")

  available_tests <- unique(gweis_res$TEST)
  message("  Available TEST values: ", paste(available_tests, collapse = ", "))

  GWAS <- gweis_res[gweis_res$TEST == int_term, ]

  if (nrow(GWAS) == 0) {
    stop(
      "No interaction term found: ", int_term, "\n",
      "Available interaction terms: ",
      paste(grep("ADDx", available_tests, value = TRUE), collapse = ", ")
    )
  }

  required_cols <- c("ID", "#CHROM", "POS", "A1", "REF", "OR", "P", "OBS_CT")
  missing_cols <- setdiff(required_cols, colnames(GWAS))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in GWEIS results: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ## ---- Convert OR to log odds for BETA ----
  GWAS$OR_numeric <- suppressWarnings(as.numeric(GWAS$OR))

  invalid_or <- is.na(GWAS$OR_numeric) | GWAS$OR_numeric <= 0
  if (any(invalid_or)) {
    n_invalid <- sum(invalid_or)
    warning("Removing ", n_invalid, " SNPs with missing or non-positive odds ratios")
    GWAS <- GWAS[!invalid_or, ]
  }

  if (nrow(GWAS) == 0) {
    stop("No valid SNPs remaining after filtering invalid odds ratios")
  }

  GWAS$BETA <- log(GWAS$OR_numeric)

  infinite_beta <- is.infinite(GWAS$BETA)
  if (any(infinite_beta)) {
    n_infinite <- sum(infinite_beta)
    warning("Removing ", n_infinite, " SNPs with infinite log odds values")
    GWAS <- GWAS[!infinite_beta, ]
  }

  if (nrow(GWAS) == 0) {
    stop("No valid SNPs remaining after filtering infinite BETA values")
  }

  ## ---- Prepare GCIM-formatted output ----
  gcim_data <- data.frame(
    order = seq_len(nrow(GWAS)),
    snpid = GWAS$ID,
    chr   = GWAS$`#CHROM`,
    bp    = GWAS$POS,
    a1    = GWAS$A1,
    a2    = GWAS$REF,
    beta  = GWAS$BETA,
    pval  = GWAS$P,
    N     = GWAS$OBS_CT,
    stringsAsFactors = FALSE
  )

  complete_cases <- complete.cases(gcim_data)
  if (!all(complete_cases)) {
    n_removed <- sum(!complete_cases)
    warning("Removing ", n_removed, " SNPs with missing values")
    gcim_data <- gcim_data[complete_cases, ]
  }

  if (nrow(gcim_data) == 0) {
    stop("No valid interaction results after filtering")
  }

  ## ---- Write GCIM-formatted file ----
  gcim_file <- file.path(tmp_dir, out_file)

  writeLines("order snpid chr bp a1 a2 beta pval N", con = gcim_file)

  write.table(
    gcim_data,
    file = gcim_file,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = " ",
    append = TRUE
  )

  ## ---- Summary output ----
  message("Binary GWEIS completed successfully")
  message("  Total SNPs tested: ", nrow(gweis_res))
  message("  Interaction effects extracted: ", nrow(gcim_data))
  message("  Interaction term: ", int_term)
  message("  Results saved to: ", tmp_dir)

  message("  Interaction effect statistics:")
  message("    Beta mean: ", round(mean(gcim_data$beta, na.rm = TRUE), 6))
  message("    Beta SD: ", round(sd(gcim_data$beta, na.rm = TRUE), 6))
  message("    Min p value: ", format(min(gcim_data$pval, na.rm = TRUE), scientific = TRUE))

  result <- list(
    glm_file = glm_file,
    gcim_file = gcim_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_interactions = nrow(gcim_data),
    interaction_term = int_term
  )

  invisible(result)
}