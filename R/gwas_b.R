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
#'   1: FID, 2: IID, 3: binary phenotype (1/2 or 0/1), 4-K: covariates.
#' @param out_prefix Character. Prefix for PLINK output inside tmp_dir (default "b_gwas").
#' @param threads Integer. Number of threads for PLINK (default 40).
#'
#' @return A list containing:
#' \describe{
#'   \item{gwas_file}{Path to PLINK .glm.logistic* file}
#'   \item{score_file}{Path to extracted additive SNP effects (ID, A1, BETA)}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{log_file}{Path to PLINK log file}
#'   \item{n_snps}{Number of additive SNP effects extracted}
#'   \item{files}{Named list of key output files (pipeline-friendly)}
#'   \item{meta}{Named list of metadata (pipeline-friendly)}
#' }
#'
#' @export
b_gwas <- function(plink_path,
                   dis_mydata,
                   dis_cov_file,
                   out_prefix = "b_gwas",
                   threads = 40) {

  # ---- Helper functions ----
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
    paste0("b_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 1 - Binary GWAS (Discovery)")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Running PLINK2 GWAS (logistic regression)...")

  # ---- Run PLINK logistic GWAS ----
  args <- c(
    "--bfile", dis_mydata,
    "--pheno", dis_cov_file,
    "--pheno-col-nums", "3",
    "--covar", dis_cov_file,
    "--covar-col-nums", "4-",
    "--glm", "hide-covar",
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

  # ---- Locate GWAS results robustly ----
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.logistic(\\..*)?$", full.names = TRUE)
  if (length(glm_files) == 0) {
    # Fallback: any glm output
    glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  }
  if (length(glm_files) == 0) {
    .stop2(
      "No PLINK GWAS output (*.glm.*) found in: ", tmp_dir, "\n",
      "Check:\n  log: ", log_file, "\n  stderr: ", stderr_file
    )
  }

  # Prefer PHENO1.glm.logistic when present
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

  # ---- Extract additive effects ----
  if (!"TEST" %in% names(gwas_res)) {
    .stop2("Column 'TEST' not found in GWAS results: ", gwas_file)
  }

  add_res <- gwas_res[gwas_res$TEST == "ADD", , drop = FALSE]
  if (nrow(add_res) == 0) {
    .stop2("No additive effects (TEST == 'ADD') found in: ", gwas_file)
  }

  if (!all(c("ID", "A1") %in% names(add_res))) {
    .stop2("Missing required columns 'ID' and/or 'A1' in GWAS results: ", gwas_file)
  }

  # ---- Derive BETA: prefer BETA column, else use log(OR) ----
  effect_source <- NULL
  n_removed <- 0
  
  if ("BETA" %in% names(add_res)) {
    effect_source <- "BETA"
    eff <- suppressWarnings(as.numeric(add_res$BETA))
    bad <- is.na(eff) | is.infinite(eff)
    n_removed <- sum(bad)
    if (any(bad)) {
      message("Removing ", n_removed, " SNPs with missing/infinite BETA values.")
      add_res <- add_res[!bad, , drop = FALSE]
      eff <- eff[!bad]
    }
    add_res$BETA <- eff
    
  } else if ("OR" %in% names(add_res)) {
    effect_source <- "log(OR)"
    orv <- suppressWarnings(as.numeric(add_res$OR))
    bad <- is.na(orv) | is.infinite(orv) | orv <= 0
    n_removed_or <- sum(bad)
    if (any(bad)) {
      message("Removing ", n_removed_or, " SNPs with missing/invalid OR values.")
      add_res <- add_res[!bad, , drop = FALSE]
      orv <- orv[!bad]
    }
    
    # Compute log(OR)
    eff <- log(orv)
    bad2 <- is.na(eff) | is.infinite(eff)
    n_removed_log <- sum(bad2)
    if (any(bad2)) {
      message("Removing ", n_removed_log, " SNPs with missing/infinite log(OR) values.")
      add_res <- add_res[!bad2, , drop = FALSE]
      eff <- eff[!bad2]
    }
    add_res$BETA <- eff
    n_removed <- n_removed_or + n_removed_log
    
  } else {
    .stop2("Neither 'BETA' nor 'OR' column found in GWAS results. File: ", gwas_file)
  }

  if (nrow(add_res) == 0) {
    .stop2("No valid SNPs remaining after filtering invalid effects. File: ", gwas_file)
  }

  # ---- Write score file for PRS construction ----
  score_file <- file.path(tmp_dir, "covariate_additive_effects.score")
  utils::write.table(
    add_res[, c("ID", "A1", "BETA"), drop = FALSE],
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
  message("Binary GWAS completed successfully!")
  message("  GWAS file:  ", gwas_file)
  message("  Score file: ", score_file)
  message("  ADD SNPs:   ", nrow(add_res))
  message("  Effect from:", effect_source)
  message("  Log file:   ", log_file)
  message(.bar())

  # ---- Return results (pipeline-friendly structure) ----
  res <- list(
    # Backwards compatible fields
    gwas_file = gwas_file,
    score_file = score_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_snps = nrow(add_res),

    # Pipeline-friendly structure for downstream steps
    files = list(
      gwas_glm = gwas_file,
      score = score_file,
      plink_log = log_file,
      stdout = stdout_file,
      stderr = stderr_file
    ),
    meta = list(
      step = "GWAS",
      trait_type = "binary",
      model = "logistic",
      pheno_col = 3L,
      covar_cols = "4-",
      threads = threads,
      test_extracted = "ADD",
      effect_source = effect_source,
      n_snps_total = nrow(gwas_res),
      n_snps_additive = nrow(add_res),
      n_snps_removed = n_removed,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}