#' Perform GWEIS for a quantitative phenotype and generate a GCIM input file
#'
#' Runs PLINK2 linear regression with SNP x covariate interaction and exports
#' interaction results in GCIM format.
#'
#' @importFrom utils read.table write.table head tail
#' @importFrom stats complete.cases
#' @param plink_path Path to PLINK2 executable.
#' @param tar_mydata Prefix for target PLINK files (.bed/.bim/.fam).
#' @param tar_pheno_file Path to phenotype file (FID IID + phenotype column).
#' @param tar_covar_file Path to covariate file (FID IID + covariates).
#' @param int_covar_index Integer index (1-based) of covariate to interact with ADD.
#' @param out_file Output GCIM file name (saved under the temporary directory).
#' @param out_prefix Prefix for PLINK output files (saved under the temporary directory).
#' @param threads Number of threads for PLINK2.
#' @param pheno_col Integer column number of phenotype in `tar_pheno_file` (>= 3).
#' @param standardize_covar Logical; standardize covariates with
#'   `--covar-variance-standardize`.
#' @param standardize_pheno Logical; standardize phenotype with
#'   `--pheno-variance-standardize`.
#' @param show_plink_output Logical; on PLINK error, print stderr head + log tail.
#' @param gcim_sep Output delimiter for GCIM file: `"space"` or `"tab"`.
#' @param keep_tmp Logical; keep the temporary working directory.
#' @return Invisibly returns a list with paths to the `.glm` file, GCIM output,
#'   temp directory, and metadata.
#' @export
q_gweis <- function(plink_path,
                    tar_mydata,
                    tar_pheno_file,
                    tar_covar_file,
                    int_covar_index = 1,
                    out_file = "gcim_prs_int.txt",
                    out_prefix = "q_gweis",
                    threads = 40,
                    pheno_col = 3,
                    standardize_covar = TRUE,
                    standardize_pheno = FALSE,
                    show_plink_output = TRUE,
                    gcim_sep = c("space", "tab"),
                    keep_tmp = TRUE)  {

  gcim_sep <- match.arg(gcim_sep)

  .bar   <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  .print_plink_debug <- function(stderr_file, log_file, n = 80) {
    if (!isTRUE(show_plink_output)) return(invisible(NULL))
    if (file.exists(stderr_file)) {
      cat("\n---- PLINK STDERR (head) ----\n")
      cat(paste(utils::head(readLines(stderr_file, warn = FALSE), n), collapse = "\n"), "\n")
    }
    if (file.exists(log_file)) {
      cat("\n---- PLINK LOG (tail) ----\n")
      cat(paste(utils::tail(readLines(log_file, warn = FALSE), n), collapse = "\n"), "\n")
    }
    invisible(NULL)
  }

  # ---- Resolve list inputs ----
  pheno_source <- NULL
  covar_source <- NULL

  if (is.list(tar_pheno_file)) {
    if (!is.null(tar_pheno_file$meta)) pheno_source <- tar_pheno_file$meta$step
    if (!is.null(tar_pheno_file$tar_pheno_file)) tar_pheno_file <- tar_pheno_file$tar_pheno_file
    else if (!is.null(tar_pheno_file$output_file)) tar_pheno_file <- tar_pheno_file$output_file
    else .stop2("tar_pheno_file list must contain tar_pheno_file/output_file")
  }

  if (is.list(tar_covar_file)) {
    if (!is.null(tar_covar_file$meta)) covar_source <- tar_covar_file$meta$step
    if (!is.null(tar_covar_file$output_file)) tar_covar_file <- tar_covar_file$output_file
    else if (!is.null(tar_covar_file$tar_covar_file)) tar_covar_file <- tar_covar_file$tar_covar_file
    else .stop2("tar_covar_file list must contain output_file/tar_covar_file")
  }

  # ---- Validate ----
  if (!file.exists(plink_path)) .stop2("PLINK not found: ", plink_path)
  for (ext in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(tar_mydata, ext))) .stop2("Missing PLINK file: ", paste0(tar_mydata, ext))
  }
  if (!file.exists(tar_pheno_file)) .stop2("Phenotype file not found: ", tar_pheno_file)
  if (!file.exists(tar_covar_file)) .stop2("Covariate file not found: ", tar_covar_file)

  int_covar_index <- as.integer(int_covar_index)
  if (is.na(int_covar_index) || int_covar_index < 1) .stop2("int_covar_index must be >= 1")

  threads <- as.integer(threads)
  if (is.na(threads) || threads < 1) threads <- 1

  pheno_col <- as.integer(pheno_col)
  if (is.na(pheno_col) || pheno_col < 3) .stop2("pheno_col must be >= 3")

  # ---- Read covar header ----
  cov_header <- tryCatch(
    utils::read.table(tar_covar_file, header = TRUE, nrows = 1,
                      stringsAsFactors = FALSE, check.names = FALSE, comment.char = ""),
    error = function(e) NULL
  )

  if (is.null(cov_header) || !all(c("FID", "IID") %in% names(cov_header))) {
    cov_header <- utils::read.table(tar_covar_file, header = FALSE, nrows = 1,
                                    stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
    if (ncol(cov_header) < 3) .stop2("Covariate file must have >=3 columns: ", tar_covar_file)
    names(cov_header) <- c("FID", "IID", paste0("COVAR", seq_len(ncol(cov_header) - 2)))
  }

  covar_names <- setdiff(names(cov_header), c("FID", "IID"))
  n_covars <- length(covar_names)
  if (n_covars < 1) .stop2("No covariates found after FID/IID in: ", tar_covar_file)
  if (int_covar_index > n_covars) .stop2("int_covar_index exceeds available covariates (", n_covars, ").")

  covar_name_used <- covar_names[int_covar_index]

  last_covar_col <- 2 + n_covars
  covar_cols <- paste(3:last_covar_col, collapse = ",")

  # ---- tmp dir ----
  tmp_dir <- file.path(tempdir(), paste0("q_gweis_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  if (!isTRUE(keep_tmp)) on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  out_pref    <- file.path(tmp_dir, out_prefix)
  log_file    <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 4 - Quantitative GWEIS (Target Sample)")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Phenotype file: ", basename(tar_pheno_file), " (pheno_col=", pheno_col, ")")
  message("Covariate file: ", basename(tar_covar_file))
  message("  Total covariates: ", n_covars)
  message("  Interaction covariate: '", covar_name_used, "' (index ", int_covar_index, ")")
  message("Covariate columns used: ", covar_cols)
  message("GCIM delimiter: ", gcim_sep)
  if (!is.null(pheno_source)) message("Phenotype source: ", pheno_source)
  if (!is.null(covar_source)) message("Covariate source: ", covar_source)

  # ---- PLINK args ----
  args <- c(
    "--bfile", tar_mydata,
    "--pheno", tar_pheno_file,
    "--pheno-col-nums", as.character(pheno_col),
    "--covar", tar_covar_file,
    "--covar-col-nums", covar_cols,
    "--glm", "interaction",
    "--allow-no-sex",
    "--threads", as.character(threads),
    "--out", out_pref
  )
  if (isTRUE(standardize_covar)) args <- c(args, "--covar-variance-standardize")
  if (isTRUE(standardize_pheno)) args <- c(args, "--pheno-variance-standardize")

  message("Running PLINK2 GWEIS (linear interaction)...")
  exit_code <- suppressWarnings(system2(plink_path, args = args, stdout = stdout_file, stderr = stderr_file))

  if (!is.null(exit_code) && exit_code != 0) {
    .print_plink_debug(stderr_file, log_file, n = 80)
    .stop2("PLINK exited with non-zero status (", exit_code, "). Check log: ", log_file)
  }

  # ---- Locate glm output ----
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.linear$", full.names = TRUE)
  if (length(glm_files) == 0) glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  if (length(glm_files) == 0) .stop2("No .glm output found. Check: ", log_file)

  glm_file <- glm_files[1]
  gweis_res <- utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE,
                                 comment.char = "", check.names = FALSE)
  if (!"TEST" %in% names(gweis_res)) .stop2("TEST column not found in: ", glm_file)

  addx_rows <- gweis_res[grepl("^ADDx", gweis_res$TEST), , drop = FALSE]
  if (nrow(addx_rows) == 0) .stop2("No ADDx* interaction rows found in: ", glm_file)

  interaction_terms <- sort(unique(addx_rows$TEST))
  preferred <- interaction_terms[interaction_terms == paste0("ADDx", covar_name_used)]
  chosen_term <- if (length(preferred) == 1) preferred else {
    if (length(interaction_terms) == 1) interaction_terms[1] else {
      if (int_covar_index > length(interaction_terms)) .stop2("int_covar_index exceeds available interaction terms.")
      interaction_terms[int_covar_index]
    }
  }

  gw <- gweis_res[gweis_res$TEST == chosen_term, , drop = FALSE]

  core_needed <- c("ID", "#CHROM", "POS", "A1", "BETA", "P", "OBS_CT")
  miss <- setdiff(core_needed, names(gw))
  if (length(miss) > 0) .stop2("Missing columns: ", paste(miss, collapse = ", "), "\nIn: ", glm_file)

  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) if (cand %in% names(gw)) { a2_col <- cand; break }
  if (is.null(a2_col)) .stop2("Cannot find allele2 column (A2/REF/ALT/ALT1) in: ", glm_file)

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
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  gcim_data <- gcim_data[stats::complete.cases(gcim_data), , drop = FALSE]
  if (nrow(gcim_data) == 0) .stop2("No valid interaction results after filtering missing values.")

  gcim_file <- file.path(tmp_dir, out_file)
  sep_out <- if (gcim_sep == "tab") "\t" else " "
  utils::write.table(gcim_data, gcim_file, quote = FALSE, col.names = TRUE,
                     row.names = FALSE, sep = sep_out)

  message(.bar())
  message("Quantitative GWEIS completed successfully!")
  message("  GLM file:  ", glm_file)
  message("  Interaction term detected: ", chosen_term)
  message("  GCIM file: ", gcim_file)
  message(.bar())

  invisible(list(
    glm_file = glm_file,
    gcim_file = gcim_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    interaction_term = chosen_term,
    covariate_name = covar_name_used,
    files = list(glm = glm_file, gcim = gcim_file, plink_log = log_file,
                 stdout = stdout_file, stderr = stderr_file),
    meta = list(step = "GWEIS", trait_type = "quantitative", int_term = chosen_term,
                int_covar_index = int_covar_index, covar_cols = covar_cols,
                gcim_sep = gcim_sep, timestamp = Sys.time())
  ))
}