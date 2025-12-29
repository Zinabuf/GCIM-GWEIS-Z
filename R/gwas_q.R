#' Perform GWAS for a quantitative covariate (discovery sample)
#'
#' Runs PLINK2 GWAS for a quantitative phenotype and extracts additive SNP effects
#' (TEST == "ADD") into a score file suitable for PRS construction.
#' All outputs are written to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @param plink_path Character. Path to PLINK2 executable.
#' @param dis_mydata Character. Prefix of PLINK binary files (.bed/.bim/.fam).
#' @param dis_cov_file Character. Phenotype + covariate file with columns:
#'   1: FID, 2: IID, 3: phenotype, 4-K: covariates.
#' @param out_prefix Character. Prefix for PLINK output inside tmp_dir (default "q_gwas").
#' @param threads Integer. Optional threads for PLINK (default 1).
#'
#' @return A list containing:
#' \describe{
#'   \item{gwas_file}{Path to PLINK .glm.linear file}
#'   \item{score_file}{Path to extracted additive SNP effects (ID, A1, BETA)}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_snps}{Number of additive SNP effects extracted}
#' }
#'
#' @export
q_gwas <- function(plink_path,
                   dis_mydata,
                   dis_cov_file,
                   out_prefix = "q_gwas",
                   threads = 1) {

  ## ---- Input validation ----
  if (!file.exists(plink_path)) {
    stop("PLINK executable not found at: ", plink_path)
  }

  for (ext in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(dis_mydata, ext))) {
      stop("Missing PLINK file: ", paste0(dis_mydata, ext))
    }
  }

  if (!file.exists(dis_cov_file)) {
    stop("Phenotype/covariate file not found: ", dis_cov_file)
  }

  ## ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("q_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)

  message("Output directory: ", tmp_dir)
  message("Step 1 of 3: Running PLINK GWAS")

  ## ---- Run PLINK GWAS (quantitative) ----
  args <- c(
    "--bfile", dis_mydata,
    "--pheno", dis_cov_file,
    "--pheno-col-nums", "3",
    "--covar", dis_cov_file,
    "--covar-col-nums", "4-n",
    "--glm",          # linear regression for quantitative phenotype
    "--allow-no-sex",
    "--covar-variance-standardize",
    "--threads", as.character(threads),
    "--out", out_pref
  )

  # Capture exit status
  exit_code <- suppressWarnings(system2(plink_path, args = args))
  log_file <- paste0(out_pref, ".log")

  if (!is.null(exit_code) && exit_code != 0) {
    warning("PLINK exited with non-zero status (", exit_code, "). Check: ", log_file)
  }

  ## ---- Locate GWAS results robustly ----
  message("Step 2 of 3: Locating GWAS results")

  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.linear$", full.names = TRUE)
  if (length(glm_files) == 0) {
    # fall back: any glm output
    glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  }
  if (length(glm_files) == 0) {
    stop("No PLINK GWAS output (*.glm.*) found in: ", tmp_dir,
         "\nCheck PLINK log at: ", log_file)
  }

  # Prefer PHENO1 if present, else take the first
  gwas_file <- glm_files[grepl("PHENO1", basename(glm_files))]
  if (length(gwas_file) == 0) gwas_file <- glm_files[1]
  gwas_file <- gwas_file[1]

  ## ---- Read GWAS results ----
  gwas_res <- tryCatch(
    utils::read.table(gwas_file,
               header = TRUE,
               stringsAsFactors = FALSE,
               comment.char = ""),
    error = function(e) {
      stop("Failed to read GWAS results from: ", gwas_file, "\n", e$message)
    }
  )

  ## ---- Extract additive SNP effects for PRS ----
  message("Step 3 of 3: Extracting additive SNP effects (TEST == 'ADD')")

  if (!"TEST" %in% names(gwas_res)) {
    stop("Column 'TEST' not found in GWAS results: ", gwas_file)
  }

  add_res <- gwas_res[gwas_res$TEST == "ADD", , drop = FALSE]
  if (nrow(add_res) == 0) {
    stop("No additive effects (TEST == 'ADD') found in: ", gwas_file)
  }

  required_cols <- c("ID", "A1", "BETA")
  missing_cols <- setdiff(required_cols, names(add_res))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in GWAS results: ",
         paste(missing_cols, collapse = ", "),
         "\nFile: ", gwas_file)
  }

  score_file <- file.path(tmp_dir, "covariate_additive_effects.score")
  utils::write.table(
    add_res[, required_cols],
    file = score_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  message("GWAS completed")
  message("  GWAS file: ", gwas_file)
  message("  Additive SNPs extracted: ", nrow(add_res))
  message("  Score file for PRS: ", score_file)

  invisible(list(
    gwas_file = gwas_file,
    score_file = score_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_snps = nrow(add_res)
  ))
}