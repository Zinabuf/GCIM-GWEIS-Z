#' Perform GWAS for a quantitative covariate (discovery sample)
#'
#' This function runs a PLINK GWAS for a quantitative phenotype and
#' extracts additive SNP effects for downstream GCIM-GWEIS-Z models.
#' All intermediate and output files are written to a temporary directory.
#'
#' @param plink_path Character. Path to the PLINK (v2.0) executable.
#' @param dis_snp Character. Prefix of PLINK binary files
#'   (.bed/.bim/.fam) for the discovery dataset.
#' @param qp_dis_cov Character. File path to phenotype and covariate file.
#'   Columns must be:
#'   \itemize{
#'     \item 1: FID
#'     \item 2: IID
#'     \item 3: Quantitative phenotype
#'     \item 4-K: Covariates
#'   }
#'
#' @return A list containing file paths to GWAS outputs stored in a
#' temporary directory:
#' \describe{
#'   \item{glm_file}{Path to PLINK .glm.linear file}
#'   \item{additive_effects}{Path to extracted additive SNP effects}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_snps}{Number of additive SNP effects extracted}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Creates a temporary directory for all outputs
#'   \item Runs PLINK GWAS with covariate adjustment
#'   \item Extracts additive SNP effects (TEST equals ADD)
#'   \item Saves results to temporary files for pipeline continuity
#' }
#'
#' @export
q_gwas <- function(plink_path,
                   dis_snp,
                   qp_dis_cov) {

  ## ---- Input validation ----
  if (!file.exists(plink_path)) {
    stop("PLINK executable not found at: ", plink_path)
  }

  if (!file.exists(paste0(dis_snp, ".bed")) ||
      !file.exists(paste0(dis_snp, ".bim")) ||
      !file.exists(paste0(dis_snp, ".fam"))) {
    stop("PLINK binary files (.bed/.bim/.fam) not found with prefix: ", dis_snp)
  }

  if (!file.exists(qp_dis_cov)) {
    stop("Phenotype and covariate file not found: ", qp_dis_cov)
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("q_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  message("Output directory: ", tmp_dir)

  out_prefix <- file.path(tmp_dir, "q_gwas")

  ## ---- Run PLINK GWAS ----
  message("Step 1 of 3: Running PLINK GWAS")

  cmd <- paste(
    shQuote(plink_path),
    "--bfile", shQuote(dis_snp),
    "--pheno", shQuote(qp_dis_cov),
    "--pheno-col-nums", 3,
    "--covar", shQuote(qp_dis_cov),
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

  glm_file <- paste0(out_prefix, ".PHENO1.glm.linear")
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

  ## ---- Extract additive SNP effects ----
  message("Step 3 of 3: Extracting additive SNP effects")

  if (!"TEST" %in% colnames(gwas_res)) {
    stop("Column 'TEST' not found in GWAS results")
  }

  add_res <- gwas_res[gwas_res$TEST == "ADD", ]

  if (nrow(add_res) == 0) {
    stop("No additive effects (TEST equals ADD) found in GWAS results")
  }

  required_cols <- c("ID", "A1", "BETA")
  missing_cols <- setdiff(required_cols, colnames(add_res))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in GWAS results: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  add_file <- file.path(tmp_dir, "covariate_additive_effects.txt")

  write.table(
    add_res[, required_cols],
    file = add_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  ## ---- Summary output ----
  message("GWAS completed successfully")
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
