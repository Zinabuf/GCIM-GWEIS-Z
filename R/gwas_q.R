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
#'   1: FID, 2: IID (or EID), 3: phenotype, 4-K: covariates (optional).
#'   Minimum 3 columns required (FID, IID/EID, phenotype).
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

  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Input validation ----
  if (!is.character(plink_path) || length(plink_path) != 1 || !nzchar(plink_path)) {
    .stop2("`plink_path` must be a non-empty character scalar.")
  }
  if (!file.exists(plink_path)) .stop2("PLINK executable not found at: ", plink_path)

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
  if (!file.exists(dis_cov_file)) .stop2("Phenotype/covariate file not found: ", dis_cov_file)

  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 1) {
    warning("Invalid `threads` value. Setting to 1.", call. = FALSE)
    threads <- 1
  }
  threads <- as.integer(threads)

  # ---- Read and validate header ----
  cov_data <- tryCatch(
    utils::read.table(
      dis_cov_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read file: ", dis_cov_file, "\n", e$message)
  )

  if (ncol(cov_data) < 3) {
    .stop2(
      "dis_cov_file must have at least 3 columns: FID, IID/EID, phenotype.\n",
      "Found only ", ncol(cov_data), " columns in: ", dis_cov_file
    )
  }

  cn <- names(cov_data)
  
  # Check for FID and IID/EID
  has_fid <- "FID" %in% cn
  has_iid <- "IID" %in% cn
  has_eid <- "EID" %in% cn
  
  if (!has_fid) {
    .stop2(
      "dis_cov_file must contain column 'FID'.\n",
      "Found: ", paste(cn, collapse = ", "), "\n",
      "File: ", dis_cov_file
    )
  }
  
  if (!has_iid && !has_eid) {
    .stop2(
      "dis_cov_file must contain column 'IID' or 'EID'.\n",
      "Found: ", paste(cn, collapse = ", "), "\n",
      "File: ", dis_cov_file
    )
  }
  
  # Rename EID to IID if needed
  if (has_eid && !has_iid) {
    message("Note: Renaming 'EID' column to 'IID' for PLINK compatibility")
    names(cov_data)[names(cov_data) == "EID"] <- "IID"
    cn <- names(cov_data)
    
    # Write temporary file with corrected names
    tmp_dir_early <- file.path(tempdir(), paste0("q_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(tmp_dir_early, recursive = TRUE, showWarnings = FALSE)
    temp_cov_file <- file.path(tmp_dir_early, "temp_cov_renamed.txt")
    utils::write.table(
      cov_data,
      file = temp_cov_file,
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t",
      quote = FALSE
    )
    dis_cov_file <- temp_cov_file
  }

  pheno_name <- cn[3]
  covar_names <- if (length(cn) > 3) cn[4:length(cn)] else character(0)
  has_covariates <- length(covar_names) > 0

  # ---- Create temp dir ----
  tmp_dir <- file.path(tempdir(), paste0("q_gwas_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 1 - Quantitative GWAS (Discovery)")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Phenotype/covariate file: ", basename(dis_cov_file))
  message("Phenotype column (3rd): ", pheno_name)
  message("Covariates: ", if (has_covariates) paste(covar_names, collapse = ", ") else "none")
  if (!has_covariates) {
    warning(
      "No covariates detected. Consider including covariates (age, sex, PCs) to control for confounding.",
      call. = FALSE
    )
  }

  # ---- Build PLINK arguments ----
  args <- c(
    "--bfile", dis_mydata,
    "--pheno", dis_cov_file,
    "--pheno-name", pheno_name
  )

  # Add covariate arguments only if covariates exist
  if (has_covariates) {
    args <- c(
      args,
      "--covar", dis_cov_file,
      "--covar-name", paste(covar_names, collapse = ","),
      "--covar-variance-standardize"
    )
  }

  # Add GLM with appropriate modifier
  glm_modifier <- if (has_covariates) "hide-covar" else "allow-no-covars"
  
  args <- c(
    args,
    "--glm", glm_modifier,
    "--allow-no-sex",
    "--threads", as.character(threads),
    "--out", out_pref
  )

  # ---- Run PLINK ----
  exit_code <- suppressWarnings(
    system2(plink_path, args = args, stdout = stdout_file, stderr = stderr_file)
  )

  if (!file.exists(log_file)) warning("Expected PLINK log not found: ", log_file, call. = FALSE)

  if (!is.null(exit_code) && exit_code != 0) {
    # Show first lines of stderr to help debugging
    stderr_head <- ""
    if (file.exists(stderr_file)) {
      x <- readLines(stderr_file, warn = FALSE)
      if (length(x) > 0) stderr_head <- paste(utils::head(x, 30), collapse = "\n")
    }
    .stop2(
      "PLINK exited with non-zero status (", exit_code, ").\n\n",
      "PLINK command:\n  ", plink_path, " ", paste(args, collapse = " "), "\n\n",
      "Check files:\n  log: ", log_file, "\n  stderr: ", stderr_file, "\n  stdout: ", stdout_file, "\n\n",
      if (nzchar(stderr_head)) paste0("stderr (first lines):\n", stderr_head, "\n") else ""
    )
  }

  # ---- Locate GWAS results ----
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.linear$", full.names = TRUE)
  if (length(glm_files) == 0) glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  if (length(glm_files) == 0) {
    .stop2(
      "No PLINK GWAS output (*.glm.*) found in: ", tmp_dir, "\n",
      "Check:\n  log: ", log_file, "\n  stderr: ", stderr_file
    )
  }

  # Look for phenotype-specific file
  pheno_pattern <- paste0(pheno_name, "\\.glm\\.linear")
  pheno_files <- glm_files[grepl(pheno_pattern, basename(glm_files), fixed = TRUE)]
  
  if (length(pheno_files) == 0) {
    # Fallback to PHENO1
    pheno_files <- glm_files[grepl("PHENO1", basename(glm_files))]
  }
  
  gwas_file <- if (length(pheno_files) > 0) pheno_files[1] else glm_files[1]
  
  if (!file.exists(gwas_file)) .stop2("Selected GWAS file does not exist: ", gwas_file)

  message("GWAS results file: ", basename(gwas_file))

  # ---- Read GWAS results ----
  gwas_res <- tryCatch(
    utils::read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) .stop2("Failed to read GWAS results from: ", gwas_file, "\n", e$message)
  )
  if (nrow(gwas_res) == 0) .stop2("GWAS file has 0 rows: ", gwas_file)
  if (!"TEST" %in% names(gwas_res)) .stop2("Column 'TEST' not found in GWAS results: ", gwas_file)

  # ---- Extract additive effects ----
  add_res <- gwas_res[gwas_res$TEST == "ADD", , drop = FALSE]
  if (nrow(add_res) == 0) .stop2("No additive effects (TEST == 'ADD') found in: ", gwas_file)

  required_cols <- c("ID", "A1", "BETA")
  missing_cols <- setdiff(required_cols, names(add_res))
  if (length(missing_cols) > 0) {
    .stop2("Missing required columns in GWAS results: ", paste(missing_cols, collapse = ", "), "\nFile: ", gwas_file)
  }

  keep <- !is.na(add_res$BETA)
  n_removed <- sum(!keep)
  if (n_removed > 0) {
    message("Removing ", n_removed, " SNPs with missing BETA values")
  }
  add_res <- add_res[keep, , drop = FALSE]
  if (nrow(add_res) == 0) .stop2("No valid SNPs remain after filtering missing BETA values.")

  # ---- Write score file ----
  score_file <- file.path(tmp_dir, "covariate_additive_effects.score")
  utils::write.table(add_res[, required_cols, drop = FALSE],
                     file = score_file, row.names = FALSE, col.names = TRUE,
                     sep = "\t", quote = FALSE)

  if (!file.exists(score_file)) .stop2("Failed to create score file: ", score_file)

  message(.bar())
  message("GWAS completed successfully!")
  message("  GWAS file:  ", gwas_file)
  message("  Score file: ", score_file)
  message("  ADD SNPs:   ", nrow(add_res))
  message("  Covariates used: ", if (has_covariates) length(covar_names) else 0)
  message("  Log file:   ", log_file)
  message(.bar())

  # ---- Return results ----
  res <- list(
    gwas_file = gwas_file,
    score_file = score_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    n_snps = nrow(add_res),
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
      pheno_name = pheno_name,
      covar_names = covar_names,
      has_covariates = has_covariates,
      n_covariates = length(covar_names),
      threads = threads,
      test_extracted = "ADD",
      n_snps_total = nrow(gwas_res),
      n_snps_additive = nrow(add_res),
      n_snps_removed = n_removed,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}