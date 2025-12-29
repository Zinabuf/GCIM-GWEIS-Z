#' Perform GWAS for a binary covariate (discovery sample)
#'
#' Runs PLINK2 logistic GWAS for a binary phenotype and extracts additive SNP
#' effects (TEST == "ADD") into a score file suitable for PRS construction.
#' Prefers PLINK's BETA column if available; otherwise uses log(OR).
#' All outputs are saved to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @param plink_path Character. Path to PLINK2 executable.
#' @param dis_mydata Character. Prefix of PLINK binary files (.bed/.bim/.fam).
#' @param dis_cov_file Character. Phenotype + covariate file with columns:
#'   1: FID, 2: IID, 3: binary phenotype, 4-K: covariates.
#' @param out_prefix Character. Prefix for PLINK output inside tmp_dir (default "b_gwas").
#' @param threads Integer. Optional threads for PLINK (default 1).
#'
#' @return A list containing:
#' \describe{
#'   \item{gwas_file}{Path to PLINK .glm.logistic* file}
#'   \item{score_file}{Path to extracted additive SNP effects (ID, A1, BETA)}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_snps}{Number of additive SNP effects extracted}
#' }
#'
#' @export
b_gwas <- function(plink_path,
                   dis_mydata,
                   dis_cov_file,
                   out_prefix = "b_gwas",
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
    paste0("b_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")

  message("Output directory: ", tmp_dir)
  message("Step 1 of 3: Running PLINK logistic GWAS")

  ## ---- Run PLINK GWAS ----
  args <- c(
    "--bfile", dis_mydata,
    "--pheno", dis_cov_file,
    "--pheno-col-nums", "3",
    "--covar", dis_cov_file,
    "--covar-col-nums", "4-n",
    "--glm", "hide-covar",  # hide covariates in output; optional
    "--allow-no-sex",
    "--covar-variance-standardize",
    "--threads", as.character(threads),
    "--out", out_pref
  )

  exit_code <- suppressWarnings(system2(plink_path, args = args))
  if (!is.null(exit_code) && exit_code != 0) {
    warning("PLINK exited with non-zero status (", exit_code, "). Check: ", log_file)
  }

  ## ---- Locate GWAS results robustly ----
  message("Step 2 of 3: Locating GWAS results")

  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.logistic(\\..*)?$", full.names = TRUE)
  if (length(glm_files) == 0) {
    # fallback to any glm output
    glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  }
  if (length(glm_files) == 0) {
    stop("No PLINK GWAS output (*.glm.*) found in: ", tmp_dir,
         "\nCheck PLINK log at: ", log_file)
  }

  # Prefer PHENO1 if present
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

  ## ---- Extract additive effects ----
  message("Step 3 of 3: Extracting additive SNP effects (TEST == 'ADD')")

  if (!"TEST" %in% names(gwas_res)) stop("Column 'TEST' not found in GWAS results: ", gwas_file)
  add_res <- gwas_res[gwas_res$TEST == "ADD", , drop = FALSE]
  if (nrow(add_res) == 0) stop("No additive effects (TEST == 'ADD') found in: ", gwas_file)

  # Determine which effect column to use
  has_beta <- "BETA" %in% names(add_res)
  has_or   <- "OR"   %in% names(add_res)

  if (!has_beta && !has_or) {
    stop("Neither 'BETA' nor 'OR' found in GWAS results. Available columns: ",
         paste(names(add_res), collapse = ", "))
  }

  if (!all(c("ID", "A1") %in% names(add_res))) {
    stop("Missing required columns 'ID' and/or 'A1' in GWAS results.")
  }

  if (has_beta) {
    beta_num <- suppressWarnings(as.numeric(add_res$BETA))
    bad <- is.na(beta_num) | is.infinite(beta_num)
    if (any(bad)) {
      warning("Removing ", sum(bad), " SNPs with missing/invalid BETA")
      add_res <- add_res[!bad, , drop = FALSE]
      beta_num <- beta_num[!bad]
    }
    add_res$BETA <- beta_num
  } else {
    or_num <- suppressWarnings(as.numeric(add_res$OR))
    bad <- is.na(or_num) | or_num <= 0
    if (any(bad)) {
      warning("Removing ", sum(bad), " SNPs with missing/non-positive OR")
      add_res <- add_res[!bad, , drop = FALSE]
      or_num <- or_num[!bad]
    }
    beta_num <- log(or_num)
    bad2 <- is.infinite(beta_num) | is.na(beta_num)
    if (any(bad2)) {
      warning("Removing ", sum(bad2), " SNPs with invalid log(OR)")
      add_res <- add_res[!bad2, , drop = FALSE]
      beta_num <- beta_num[!bad2]
    }
    add_res$BETA <- beta_num
  }

  if (nrow(add_res) == 0) stop("No valid SNPs remaining after filtering invalid effects.")

  ## ---- Save score file for PRS ----
  score_file <- file.path(tmp_dir, "covariate_additive_effects.score")
  utils::write.table(
    add_res[, c("ID", "A1", "BETA")],
    file = score_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  message("Logistic GWAS completed")
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