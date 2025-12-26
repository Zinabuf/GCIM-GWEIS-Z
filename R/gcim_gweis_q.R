#' Perform GWEIS for a quantitative phenotype and generate a GCIM input file
#'
#' Runs PLINK2 linear regression with SNP x covariate interaction and exports
#' interaction results (ADDxCOVARk) in GCIM format.
#'
#' @importFrom utils read.table write.table
#' @importFrom stats complete.cases
#' @param plink_path Character. Path to PLINK2 executable.
#' @param tar_mydata Character. Prefix of PLINK binary files (.bed/.bim/.fam).
#' @param tar_pheno_file Character. Phenotype file with columns: FID IID PHENO
#'   (PHENO should be in column 3).
#' @param tar_covar_file Character. Covariate file with columns: FID IID COVAR1...
#'   (COVAR1 is column 3, used as COVAR1 in PLINK naming).
#' @param int_covar_index Integer. Covariate index for interaction term:
#'   1 -> ADDxCOVAR1, 2 -> ADDxCOVAR2, ...
#' @param out_file Character. Name of GCIM interaction output file.
#' @param out_prefix Character. PLINK output prefix inside tmp_dir (default "q_gweis").
#' @param threads Integer. Optional threads for PLINK (default 40).
#'
#' @return A list with file paths and metadata.
#' @export
q_gweis <- function(plink_path,
                    tar_mydata,
                    tar_pheno_file,
                    tar_covar_file,
                    int_covar_index = 1,
                    out_file = "gcim_prs_int.txt",
                    out_prefix = "q_gweis",
                    threads = 40) {

  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Input validation ----
  if (!file.exists(plink_path)) .stop2("PLINK executable not found at: ", plink_path)

  for (ext in c(".bed", ".bim", ".fam")) {
    fp <- paste0(tar_mydata, ext)
    if (!file.exists(fp)) .stop2("Missing PLINK file: ", fp)
  }

  if (!file.exists(tar_pheno_file)) .stop2("Phenotype file not found: ", tar_pheno_file)
  if (!file.exists(tar_covar_file)) .stop2("Covariate file not found: ", tar_covar_file)

  if (!is.numeric(int_covar_index) || length(int_covar_index) != 1 || is.na(int_covar_index) || int_covar_index < 1) {
    .stop2("int_covar_index must be a positive integer")
  }
  int_covar_index <- as.integer(int_covar_index)

  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 1) {
    warning("Invalid `threads` value. Setting to 1.", call. = FALSE)
    threads <- 1
  }
  threads <- as.integer(threads)

  # ---- Check covariate header / count + name of selected covariate ----
  cov_header <- tryCatch(
    utils::read.table(tar_covar_file, header = TRUE, nrows = 1, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) .stop2("Failed to read covariate file header: ", tar_covar_file, "\n", e$message)
  )

  if (!all(c("FID", "IID") %in% names(cov_header))) {
    .stop2("Covariate file must contain 'FID' and 'IID' columns")
  }

  covar_names <- setdiff(names(cov_header), c("FID", "IID"))
  n_covars <- length(covar_names)
  if (n_covars < 1) .stop2("Covariate file must contain at least one covariate after FID and IID")
  if (int_covar_index > n_covars) {
    .stop2("int_covar_index (", int_covar_index, ") exceeds number of covariates (", n_covars, ")")
  }
  covar_name_used <- covar_names[int_covar_index]

  # ---- Temp dir ----
  tmp_dir <- file.path(tempdir(), paste0("q_gweis_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_pref <- file.path(tmp_dir, out_prefix)
  log_file <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  int_term <- paste0("ADDxCOVAR", int_covar_index)

  message(.bar())
  message("GCIM-GWEIS-Z: Quantitative GWEIS (linear interaction)")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Interaction term: ", int_term, "  (covariate: ", covar_name_used, ")")

  # ---- Run PLINK2 GWEIS ----
  args <- c(
    "--bfile", tar_mydata,
    "--pheno", tar_pheno_file,
    "--pheno-col-nums", "3",
    "--covar", tar_covar_file,
    "--covar-col-nums", "3-n",
    "--glm", "interaction", "hide-covar",
    "--allow-no-sex",
    "--covar-variance-standardize",
    "--pheno-variance-standardize",
    "--threads", as.character(threads),
    "--out", out_pref
  )

  exit_code <- suppressWarnings(system2(plink_path, args = args, stdout = stdout_file, stderr = stderr_file))
  if (!is.null(exit_code) && exit_code != 0) {
    .stop2(
      "PLINK exited with non-zero status (", exit_code, ").\n",
      "Check:\n  log: ", log_file, "\n  stderr: ", stderr_file, "\n  stdout: ", stdout_file
    )
  }

  # ---- Locate GWEIS results ----
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.linear$", full.names = TRUE)
  if (length(glm_files) == 0) {
    .stop2("No .glm.linear output found in: ", tmp_dir, "\nCheck PLINK log at: ", log_file)
  }
  pheno1 <- glm_files[grepl("PHENO1", basename(glm_files))]
  glm_file <- if (length(pheno1) > 0) pheno1[1] else glm_files[1]

  # ---- Read GWEIS results ----
  gweis_res <- tryCatch(
    utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE),
    error = function(e) .stop2("Failed to read GWEIS results from ", glm_file, "\n", e$message)
  )

  if (!"TEST" %in% names(gweis_res)) .stop2("Column 'TEST' not found in GWEIS output: ", glm_file)

  # ---- Extract interaction effects ----
  gw <- gweis_res[gweis_res$TEST == int_term, , drop = FALSE]
  if (nrow(gw) == 0) {
    tests <- unique(gweis_res$TEST)
    .stop2(
      "No rows found for interaction term: ", int_term, "\n",
      "Available interaction terms: ",
      paste(grep("^ADDxCOVAR", tests, value = TRUE), collapse = ", ")
    )
  }

  # ---- Required columns ----
  core_needed <- c("ID", "#CHROM", "POS", "A1", "BETA", "P", "OBS_CT")
  missing_core <- setdiff(core_needed, names(gw))
  if (length(missing_core) > 0) {
    .stop2("Missing required columns in GWEIS results: ",
           paste(missing_core, collapse = ", "),
           "\nColumns available: ", paste(names(gw), collapse = ", "))
  }

  # Allele 2 column
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(gw)) { a2_col <- cand; break }
  }
  if (is.null(a2_col)) {
    .stop2("Could not find a second allele column (tried A2/REF/ALT/ALT1). Columns: ",
           paste(names(gw), collapse = ", "))
  }

  # ---- Build GCIM data ----
  gcim_data <- data.frame(
    order = seq_len(nrow(gw)),
    snpid = gw$ID,
    chr   = gw$`#CHROM`,
    bp    = gw$POS,
    a1    = gw$A1,
    a2    = gw[[a2_col]],
    beta  = gw$BETA,
    pval  = gw$P,
    N     = gw$OBS_CT,
    stringsAsFactors = FALSE
  )

  keep <- stats::complete.cases(gcim_data)
  if (!all(keep)) {
    warning("Removing ", sum(!keep), " SNPs with missing values", call. = FALSE)
    gcim_data <- gcim_data[keep, , drop = FALSE]
  }
  if (nrow(gcim_data) == 0) .stop2("No valid interaction results after filtering missing values.")

  # ---- Write GCIM file ----
  gcim_file <- file.path(tmp_dir, out_file)
  writeLines("order snpid chr bp a1 a2 beta pval N", con = gcim_file)
  utils::write.table(gcim_data, file = gcim_file, quote = FALSE, col.names = FALSE,
                     row.names = FALSE, sep = " ", append = TRUE)

  message("GWEIS completed")
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
    interaction_term = int_term,
    covariate_name = covar_name_used,
    files = list(
      glm = glm_file,
      gcim = gcim_file,
      plink_log = log_file,
      stdout = stdout_file,
      stderr = stderr_file
    ),
    meta = list(
      step = "GWEIS",
      trait_type = "quantitative",
      model = "linear_interaction",
      int_covar_index = int_covar_index,
      int_term = int_term,
      covariate_name = covar_name_used,
      threads = threads,
      a2_col = a2_col
    )
  ))
}
