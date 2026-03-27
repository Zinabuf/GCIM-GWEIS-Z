#' Perform GWEIS for a binary phenotype and generate a GCIM input file
#'
#' Runs PLINK2 logistic/Firth regression with SNP x covariate interaction and
#' exports the selected interaction result in GCIM format.
#'
#' Key behavior:
#' - Keeps all covariates in the model as adjustment covariates
#' - Includes only ONE SNP x covariate interaction term via --parameters
#' - Writes GCIM file as tab-delimited by default
#' - Retains richer GCIM output (chr, bp, alleles, beta, p, N)
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
#' @param standardize_covar Logical; standardize covariates using
#'   `--covar-variance-standardize`.
#' @param ensure_12 Logical; if TRUE, recode binary phenotype 0/1 -> 1/2 when needed.
#' @param firth_mode One of "default", "firth-fallback", "firth", or "no-firth".
#' @param show_plink_output Logical; on PLINK error, print stderr head + log tail.
#' @param gcim_sep Output delimiter for GCIM file: `"tab"` or `"space"`.
#' @param keep_tmp Logical; keep the temporary working directory.
#' @return Invisibly returns a list containing paths to the `.glm` file, GCIM file,
#'   temporary directory, and metadata.
#'
#' @export
b_gweis <- function(plink_path,
                    tar_mydata,
                    tar_pheno_file,
                    tar_covar_file,
                    int_covar_index = 1,
                    out_file = "gcim_prs_int.txt",
                    out_prefix = "b_gweis",
                    threads = 40,
                    pheno_col = 3,
                    standardize_covar = TRUE,
                    ensure_12 = TRUE,
                    firth_mode = c("default", "firth-fallback", "firth", "no-firth"),
                    show_plink_output = TRUE,
                    gcim_sep = c("tab", "space"),
                    keep_tmp = TRUE) {

  gcim_sep   <- match.arg(gcim_sep)
  firth_mode <- match.arg(firth_mode)

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

  .read_table_flex <- function(path, min_cols = 3) {
    x <- tryCatch(
      utils::read.table(path, header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE, comment.char = ""),
      error = function(e) NULL
    )

    ok <- !is.null(x) && ncol(x) >= 2 && any(c("FID", "IID", "EID") %in% names(x))
    if (!ok) {
      x <- utils::read.table(path, header = FALSE, stringsAsFactors = FALSE,
                             check.names = FALSE, comment.char = "")
      if (nrow(x) == 0) .stop2("File has 0 rows: ", path)
      if (ncol(x) < min_cols) .stop2("File must have at least ", min_cols, " columns: ", path)
      names(x)[1:2] <- c("FID", "IID")
    }

    if ("EID" %in% names(x) && !("IID" %in% names(x))) {
      names(x)[names(x) == "EID"] <- "IID"
    }

    if (!all(c("FID", "IID") %in% names(x))) {
      .stop2("File must contain FID and IID columns: ", path,
             "\nColumns: ", paste(names(x), collapse = ", "))
    }

    if (ncol(x) < min_cols) .stop2("File must have at least ", min_cols, " columns: ", path)
    x
  }

  .validate_binary_and_maybe_recode <- function(pheno_df, pheno_col, tmp_dir, ensure_12) {
    v <- suppressWarnings(as.numeric(pheno_df[[pheno_col]]))
    if (all(is.na(v))) .stop2("Phenotype column ", pheno_col, " is entirely NA/non-numeric.")

    u <- sort(unique(v[!is.na(v)]))
    if (length(u) < 2) .stop2("Binary phenotype has <2 unique values: ", paste(u, collapse = ", "))
    if (length(u) > 2) .stop2("Phenotype column appears non-binary. Unique values: ",
                              paste(utils::head(u, 10), collapse = ", "))

    if (isTRUE(ensure_12) && identical(u, c(0, 1))) {
      ph2 <- pheno_df
      ph2[[pheno_col]] <- ifelse(is.na(v), NA, ifelse(v == 0, 1, 2))
      out <- file.path(tmp_dir, "pheno_recode_12.tsv")
      utils::write.table(ph2, out, row.names = FALSE, col.names = TRUE,
                         sep = "\t", quote = FALSE)
      return(list(path = out, recoded = TRUE, levels = c(1, 2)))
    }

    list(path = NULL, recoded = FALSE, levels = u)
  }

  # ---- Resolve list inputs (pipeline output) ----
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

  # ---- Read covariate header / covariates ----
  cov_df <- .read_table_flex(tar_covar_file, min_cols = 3)
  covar_names <- setdiff(names(cov_df), c("FID", "IID"))
  n_covars <- length(covar_names)

  if (n_covars < 1) .stop2("No covariates found after FID/IID in: ", tar_covar_file)
  if (int_covar_index > n_covars) .stop2("int_covar_index exceeds available covariates (", n_covars, ").")

  covar_name_used <- covar_names[int_covar_index]
  covar_cols <- paste(3:(n_covars + 2), collapse = ",")

  # ---- Build parameter selection from Function B logic ----
  # Predictor order under --glm interaction:
  # 1 = ADD
  # 2:(n_covars+1) = covariate main effects
  # (n_covars+2):(2*n_covars+1) = ADD x covariate interactions
  add_idx <- 1
  covar_idx <- 2:(n_covars + 1)
  interaction_idx <- (n_covars + 1) + int_covar_index

  parameters_keep <- c(add_idx, covar_idx, interaction_idx)
  parameters_keep_str <- paste(parameters_keep, collapse = ",")

  # ---- tmp dir ----
  tmp_dir <- file.path(tempdir(), paste0("b_gweis_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  if (!isTRUE(keep_tmp)) on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  out_pref    <- file.path(tmp_dir, out_prefix)
  log_file    <- paste0(out_pref, ".log")
  stdout_file <- paste0(out_pref, ".stdout")
  stderr_file <- paste0(out_pref, ".stderr")

  # ---- Read phenotype robustly + validate binary + optional recode ----
  pheno_df <- .read_table_flex(tar_pheno_file, min_cols = pheno_col)
  rec <- .validate_binary_and_maybe_recode(pheno_df, pheno_col, tmp_dir, ensure_12)
  pheno_path_used <- if (isTRUE(rec$recoded)) rec$path else tar_pheno_file

  message(.bar())
  message("GCIM-GWEIS-Z: Step 4 - Binary GWEIS (Target Sample)")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Phenotype file: ", basename(pheno_path_used), " (pheno_col=", pheno_col, ")")
  message("  Binary levels observed: ", paste(rec$levels, collapse = ", "),
          if (rec$recoded) " (recoded 0/1 -> 1/2)" else "")
  message("Covariate file: ", basename(tar_covar_file))
  message("  Total covariates: ", n_covars)
  message("  Interaction covariate: '", covar_name_used, "' (index ", int_covar_index, ")")
  message("Covariate columns used: ", covar_cols)
  message("PLINK --parameters keep: ", parameters_keep_str)
  message("GCIM delimiter: ", gcim_sep)
  message("Firth mode: ", firth_mode)
  if (!is.null(pheno_source)) message("Phenotype source: ", pheno_source)
  if (!is.null(covar_source)) message("Covariate source: ", covar_source)

  # ---- PLINK args ----
  args <- c(
    "--bfile", tar_mydata,
    "--pheno", pheno_path_used,
    "--pheno-col-nums", as.character(pheno_col),
    "--covar", tar_covar_file,
    "--covar-col-nums", covar_cols,
    "--glm", "interaction",
    "--parameters", parameters_keep_str,
    "--allow-no-sex"
  )

  if (isTRUE(standardize_covar)) {
    args <- c(args, "--covar-variance-standardize")
  }

  if (firth_mode != "default") {
    args <- c(args, paste0("--", firth_mode))
  }

  args <- c(args,
            "--threads", as.character(threads),
            "--out", out_pref)

  message("Running PLINK2 GWEIS (logistic/Firth interaction)...")
  exit_code <- suppressWarnings(system2(plink_path, args = args,
                                        stdout = stdout_file, stderr = stderr_file))

  if (!is.null(exit_code) && exit_code != 0) {
    .print_plink_debug(stderr_file, log_file, n = 80)
    .stop2("PLINK exited with non-zero status (", exit_code, "). Check log: ", log_file)
  }

  # ---- Locate glm output ----
  glm_files <- list.files(tmp_dir, pattern = "\\.glm\\.", full.names = TRUE)
  if (length(glm_files) == 0) .stop2("No .glm output found. Check: ", log_file)

  ph1 <- glm_files[grepl("PHENO1", basename(glm_files))]
  glm_file <- if (length(ph1) > 0) ph1[1] else glm_files[1]

  gweis_res <- utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE,
                                 comment.char = "", check.names = FALSE)
  if (!"TEST" %in% names(gweis_res)) .stop2("TEST column not found in: ", glm_file)

  # ---- Select interaction term robustly ----
  interaction_rows <- gweis_res[grepl("^ADDx", gweis_res$TEST), , drop = FALSE]
  if (nrow(interaction_rows) == 0) {
    .stop2("No interaction rows found in PLINK output.\n",
           "This may indicate parameter indexing mismatch or PLINK naming differences.")
  }

  preferred_patterns <- c(
    paste0("^ADDx", covar_name_used, "$"),
    paste0("^ADDxCOVAR", int_covar_index, "$"),
    paste0("^ADDx", int_covar_index, "$")
  )

  chosen_term <- NULL
  for (pat in preferred_patterns) {
    hit <- interaction_rows$TEST[grepl(pat, interaction_rows$TEST)]
    if (length(hit) >= 1) {
      chosen_term <- hit[1]
      break
    }
  }

  if (is.null(chosen_term)) {
    if (length(unique(interaction_rows$TEST)) == 1) {
      chosen_term <- unique(interaction_rows$TEST)[1]
    } else {
      .stop2("Could not uniquely identify the intended interaction term.\n",
             "Observed interaction TEST values: ",
             paste(unique(interaction_rows$TEST), collapse = ", "))
    }
  }

  gw <- interaction_rows[interaction_rows$TEST == chosen_term, , drop = FALSE]
  if (nrow(gw) == 0) .stop2("Chosen interaction term not found after selection: ", chosen_term)

  # ---- Required columns (keep Function A richer output) ----
  core_needed <- c("ID", "#CHROM", "POS", "A1", "P", "OBS_CT")
  miss <- setdiff(core_needed, names(gw))
  if (length(miss) > 0) .stop2("Missing columns: ", paste(miss, collapse = ", "), "\nIn: ", glm_file)

  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(gw)) {
      a2_col <- cand
      break
    }
  }
  if (is.null(a2_col)) .stop2("Cannot find allele2 column (A2/REF/ALT/ALT1) in: ", glm_file)

  # ---- Effect size: keep Function A fallback logic ----
  effect_source <- NULL
  if ("BETA" %in% names(gw)) {
    effect_source <- "BETA"
    beta <- suppressWarnings(as.numeric(gw$BETA))
  } else if ("OR" %in% names(gw)) {
    effect_source <- "log(OR)"
    beta <- log(suppressWarnings(as.numeric(gw$OR)))
  } else if ("OR_ALT" %in% names(gw)) {
    effect_source <- "log(OR_ALT)"
    beta <- log(suppressWarnings(as.numeric(gw$OR_ALT)))
  } else {
    .stop2("No effect size column found (expected BETA, OR, or OR_ALT).")
  }

  if (any(!is.finite(beta))) {
    bad <- !is.finite(beta)
    warning("Dropping ", sum(bad), " SNPs with invalid effect values (", effect_source, ").",
            call. = FALSE)
    gw <- gw[!bad, , drop = FALSE]
    beta <- beta[!bad]
  }
  if (nrow(gw) == 0) .stop2("No valid interaction results remain after filtering invalid effects.")

  gcim_data <- data.frame(
    order = seq_len(nrow(gw)),
    snpid = gw$ID,
    chr   = gw$`#CHROM`,
    bp    = gw$POS,
    a1    = gw$A1,
    a2    = gw[[a2_col]],
    beta  = beta,
    pval  = gw$P,
    N     = gw$OBS_CT,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  gcim_data <- gcim_data[stats::complete.cases(gcim_data), , drop = FALSE]
  if (nrow(gcim_data) == 0) .stop2("No valid interaction results after filtering missing values.")

  gcim_file <- file.path(tmp_dir, out_file)
  sep_out <- if (gcim_sep == "tab") "\t" else " "
  utils::write.table(gcim_data, file = gcim_file, quote = FALSE,
                     col.names = TRUE, row.names = FALSE, sep = sep_out)

  message(.bar())
  message("Binary GWEIS completed successfully!")
  message("  GLM file:  ", glm_file)
  message("  Interaction term detected: ", chosen_term)
  message("  Effect measure: ", effect_source)
  message("  GCIM file: ", gcim_file)
  message(.bar())

  invisible(list(
    glm_file = glm_file,
    gcim_file = gcim_file,
    tmp_dir = tmp_dir,
    log_file = log_file,
    interaction_term = chosen_term,
    covariate_name = covar_name_used,
    parameters_keep = parameters_keep,
    interaction_idx = interaction_idx,
    files = list(
      glm = glm_file,
      gcim = gcim_file,
      plink_log = log_file,
      stdout = stdout_file,
      stderr = stderr_file,
      pheno_used = pheno_path_used
    ),
    meta = list(
      step = "GWEIS",
      trait_type = "binary",
      model = "logistic_or_firth_interaction",
      int_covar_index = int_covar_index,
      int_term = chosen_term,
      interaction_idx = interaction_idx,
      parameters_keep = parameters_keep,
      covar_cols = covar_cols,
      gcim_sep = gcim_sep,
      firth_mode = firth_mode,
      effect_source = effect_source,
      timestamp = Sys.time()
    )
  ))
}