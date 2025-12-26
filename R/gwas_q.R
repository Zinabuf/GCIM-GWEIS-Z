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
#' @param threads Integer. Number of threads for PLINK (default 40).
#'
#' @return A list containing:
#' \describe{
#'   \item{gwas_file}{Path to PLINK .glm.linear file}
#'   \item{score_file}{Path to extracted additive SNP effects (ID, A1, BETA)}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_snps}{Number of additive SNP effects extracted}
#'   \item{files}{Named list of key output files (pipeline-friendly)}
#'   \item{meta}{Named list of metadata (pipeline-friendly)}
#' }
#'
#' @export
q_gwas <- function(plink_path,
                   dis_mydata,
                   dis_cov_file,
                   out_prefix = "q_gwas",
                   threads = 40) {

  # ---- helpers ----
  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Input validation ----
  if (!is.character(plink_path) || length(plink_path) != 1 || !nzchar(plink_path)) {
    .stop2("`plink_path` must be a non-empty character scalar.")
  }
  if (!file.exists(plink_path)) {
    .stop2("PLINK executable not found at: ", plink_path)
  }

  if (!is.character(dis_mydata) || length(dis_mydata) != 1 || !nzchar(dis_mydata)) {
    .stop2("`dis_mydata` must be a non-empty character scalar prefix to .bed/.bim/.fam.")
  }
  for (ext in c(".bed", ".bim", ".fam")) {
    fp <- paste0(dis_mydata, ext)
    if (!file.exists(fp)) .stop2("Missing PLINK file: ", fp)
  }

  if (!is.character(dis_cov_file) || length(dis_cov_file) != 1 || !nzchar(dis_cov_file)) {
    .stop2("`dis_cov_file` must be a non-empty character scalar.")
  }
  if (!file.exists(dis_cov_file)) {
    .stop2("Phenotype/covariate file not found: ", dis_cov_file)
  }

  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 1) {
    warning("Invalid `threads` value. Setting to 1.", call. = FALSE)
    threads <- 1
  }
  threads <- as.integer(threads)

  # ---- Create step-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("q_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 1 - Quantitative GWAS")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Running PLINK2 GWAS (linear)...")

  # ---- Run PLINK GWAS (quantitative) ----
  args <- c(
    "--bfile", dis_mydata,
    "--pheno", dis_cov_file,
    "--pheno-col-nums", "3",
    "--covar", dis_cov_file,
    "--covar-col-nums", "4-n",
    "--glm", "hide-covar",
    "--allow-no-sex",
    "--covar-variance-standardize",
    "--threads", as.character(threads),
    "--out", out_pref
  )

  exit_code <- suppressWarnings(
    system2(plink_path, args = args, stdout = stdout_file, stderr = stderr_file)
  )

  # If PLINK didn't create a log, that's usually a strong signal of failure
  if (!file.exists(log_file)) {
    warning("Expected PLINK log not found: ", log_file, call. = FALSE)
  }

  if (!is.null(exit_code) && exit_code != 0) {
    .stop2(
      "PLINK exited with non-zero status (", exit_code, ").\n",
      "Check:\n  log: ", log_file, "\n  stderr: ", stderr_file, "\n  stdout: ", stdout_file
    )
  }

  # ---- Locate GWAS results robustly ----
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.linear$", full.names = TRUE)
  if (length(glm_files) == 0) {
    # fallback: any glm output if linear suffix differs in a future build
    glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  }
  if (length(glm_files) == 0) {
    .stop2(
      "No PLINK GWAS output (*.glm.*) found in: ", tmp_dir, "\n",
      "Check:\n  log: ", log_file, "\n  stderr: ", stderr_file
    )
  }

  # Prefer PHENO1.glm.linear when present
  pheno1 <- glm_files[grepl("PHENO1", basename(glm_files))]
  gwas_file <- if (length(pheno1) > 0) pheno1[1] else glm_files[1]

  if (!file.exists(gwas_file)) {
    .stop2("Selected GWAS file does not exist: ", gwas_file)
  }

  message("GWAS results file: ", basename(gwas_file))

  # ---- Read GWAS results ----
  gwas_res <- tryCatch(
    utils::read.table(
      gwas_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read GWAS results from: ", gwas_file, "\n", e$message)
  )

  if (nrow(gwas_res) == 0) {
    .stop2("GWAS file has 0 rows: ", gwas_file)
  }

  # ---- Extract additive SNP effects for PRS ----
  if (!"TEST" %in% names(gwas_res)) {
    .stop2("Column 'TEST' not found in GWAS results: ", gwas_file)
  }

  add_res <- gwas_res[gwas_res$TEST == "ADD", , drop = FALSE]
  if (nrow(add_res) == 0) {
    .stop2("No additive effects (TEST == 'ADD') found in: ", gwas_file)
  }

  required_cols <- c("ID", "A1", "BETA")
  missing_cols <- setdiff(required_cols, names(add_res))
  if (length(missing_cols) > 0) {
    .stop2(
      "Missing required columns in GWAS results: ",
      paste(missing_cols, collapse = ", "),
      "\nFile: ", gwas_file
    )
  }

  # Drop missing BETA
  keep <- !is.na(add_res$BETA)
  if (!all(keep)) {
    message("Removing ", sum(!keep), " SNPs with missing BETA.")
    add_res <- add_res[keep, , drop = FALSE]
  }

  score_file <- file.path(tmp_dir, "covariate_additive_effects.score")
  utils::write.table(
    add_res[, required_cols, drop = FALSE],
    file = score_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  if (!file.exists(score_file)) {
    .stop2("Failed to create score file: ", score_file)
  }

  message(.bar())
  message("GWAS completed successfully!")
  message("  GWAS file:  ", gwas_file)
  message("  Score file: ", score_file)
  message("  ADD SNPs:   ", nrow(add_res))
  message("  Log file:   ", log_file)
  message(.bar())

  # ---- return (pipeline-friendly + backwards compatible fields) ----
  res <- list(
    # backwards compatible fields
    gwas_file = gwas_file,
    score_file = score_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_snps = nrow(add_res),

    # pipeline-friendly structure for downstream steps
    files = list(
      gwas_glm = gwas_file,
      score = score_file,
      plink_log = log_file,
      stdout = stdout_file,
      stderr = stderr_file
    ),
    meta = list(
      step = "GWAS",
      trait_type = "quantitative",
      pheno_col = 3L,
      covar_cols = "4-n",
      threads = threads,
      test_extracted = "ADD"
    )
  )

  invisible(res)
}