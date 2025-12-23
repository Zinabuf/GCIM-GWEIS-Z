#' Perform GWAS for a binary covariate (discovery sample)
#'
#' This function runs a PLINK logistic GWAS for a binary phenotype and
#' extracts additive SNP effects (log odds) for use in downstream
#' GCIM-GWEIS-Z models. All outputs are saved to a temporary directory.
#'
#' @param plink_path Character. Path to the PLINK (v2.0) executable.
#' @param dis_snp Character. Prefix of PLINK binary files
#'   (.bed/.bim/.fam) for the discovery dataset.
#' @param bp_dis_cov Character. File path to phenotype and covariate file.
#'   Columns must be:
#'   \itemize{
#'     \item 1: FID
#'     \item 2: IID
#'     \item 3: Binary phenotype (0/1 or 1/2)
#'     \item 4-K: Covariates
#'   }
#'
#' @return A list containing file paths to GWAS outputs stored in a
#' temporary directory:
#' \describe{
#'   \item{glm_file}{Path to PLINK .glm.logistic.hybrid file}
#'   \item{additive_effects}{Path to extracted additive SNP log odds}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_snps}{Number of SNPs with additive effects}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Creates a temporary directory for all outputs
#'   \item Runs PLINK logistic GWAS with covariate adjustment
#'   \item Extracts additive SNP effects (TEST equals ADD)
#'   \item Converts odds ratios to log odds for BETA values
#'   \item Saves results to temporary files for pipeline continuity
#' }
#'
#' @export
b_gwas <- function(plink_path,
                   dis_snp,
                   bp_dis_cov) {

  ## ---- Input validation ----
  if (!file.exists(plink_path)) {
    stop("PLINK executable not found at: ", plink_path)
  }

  if (!file.exists(paste0(dis_snp, ".bed")) ||
      !file.exists(paste0(dis_snp, ".bim")) ||
      !file.exists(paste0(dis_snp, ".fam"))) {
    stop("PLINK binary files (.bed/.bim/.fam) not found with prefix: ", dis_snp)
  }

  if (!file.exists(bp_dis_cov)) {
    stop("Phenotype and covariate file not found: ", bp_dis_cov)
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("b_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  message("Output directory: ", tmp_dir)

  out_prefix <- file.path(tmp_dir, "b_gwas")

  ## ---- Run PLINK logistic GWAS ----
  message("Step 1 of 3: Running PLINK logistic GWAS")

  cmd <- paste(
    shQuote(plink_path),
    "--bfile", shQuote(dis_snp),
    "--pheno", shQuote(bp_dis_cov),
    "--pheno-col-nums", 3,
    "--covar", shQuote(bp_dis_cov),
    "--covar-col-nums", "4-n",
    "--glm",
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

  ## ---- Locate GWAS results ----
  message("Step 2 of 3: Locating GWAS results")

  glm_file <- paste0(out_prefix, ".PHENO1.glm.logistic.hybrid")
  log_file <- paste0(out_prefix, ".log")

  if (!file.exists(glm_file)) {
    stop(
      "GWAS output file not found: ", glm_file,
      "\nCheck PLINK log at: ", log_file
    )
  }

  ## ---- Read GWAS results ----
  gwas_res <- tryCatch(
    {
      read.table(
        glm_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        comment.char = ""
      )
    },
    error = function(e) {
      stop("Failed to read GWAS results from ", glm_file, "\n", e$message)
    }
  )

  ## ---- Extract additive effects ----
  message("Step 3 of 3: Extracting additive SNP effects")

  if (!"TEST" %in% colnames(gwas_res)) {
    stop("Column 'TEST' not found in GWAS results")
  }

  add_res <- gwas_res[gwas_res$TEST == "ADD", ]

  if (nrow(add_res) == 0) {
    stop("No additive effects (TEST equals ADD) found in GWAS results")
  }

  required_cols <- c("ID", "A1", "OR")
  missing_cols <- setdiff(required_cols, colnames(add_res))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in GWAS results: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ## ---- Convert OR to log odds ----
  add_res$OR_numeric <- suppressWarnings(as.numeric(add_res$OR))

  invalid_or <- is.na(add_res$OR_numeric) | add_res$OR_numeric <= 0
  if (any(invalid_or)) {
    n_invalid <- sum(invalid_or)
    warning(
      "Removing ", n_invalid,
      " SNPs with missing or non-positive odds ratios"
    )
    add_res <- add_res[!invalid_or, ]
  }

  if (nrow(add_res) == 0) {
    stop("No valid SNPs remaining after filtering invalid odds ratios")
  }

  add_res$BETA <- log(add_res$OR_numeric)

  infinite_beta <- is.infinite(add_res$BETA)
  if (any(infinite_beta)) {
    n_infinite <- sum(infinite_beta)
    warning(
      "Removing ", n_infinite,
      " SNPs with infinite log odds values"
    )
    add_res <- add_res[!infinite_beta, ]
  }

  if (nrow(add_res) == 0) {
    stop("No valid SNPs remaining after filtering infinite BETA values")
  }

  ## ---- Save additive effects ----
  add_file <- file.path(tmp_dir, "covariate_additive_effects.txt")

  write.table(
    add_res[, c("ID", "A1", "BETA")],
    file = add_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  ## ---- Summary output ----
  message("Logistic GWAS completed successfully")
  message("  Total SNPs tested: ", nrow(gwas_res))
  message("  Additive effects extracted: ", nrow(add_res))
  message("  Results saved to: ", tmp_dir)

  ## ---- Return paths for pipeline continuity ----
  result <- list(
    glm_file = glm_file,
    additive_effects = add_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_snps = nrow(add_res)
  )

  invisible(result)
}