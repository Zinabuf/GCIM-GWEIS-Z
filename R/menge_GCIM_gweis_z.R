#' Munge interaction summary statistics for LDSC
#'
#' Takes a PLINK2 `.glm` file (linear or logistic), extracts rows for a specified
#' interaction term (e.g., "ADDxCOVAR1"), creates an LDSC-ready table
#' (SNP, A1, A2, N, Z), and munges it with `munge_sumstats.py`.
#' All outputs are saved to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @importFrom stats complete.cases
#' @param munge_path Character. Path to `munge_sumstats.py`.
#' @param glm_file Character OR list. Path to PLINK `.glm.*` file OR list from
#'   `q_gweis()`/`b_gweis()` containing `$glm_file`.
#' @param interaction_term Character OR NULL. Interaction test name, e.g. `"ADDxCOVAR1"`.
#'   If NULL and a list is provided, attempts to extract from metadata.
#' @param hm3_snplist Character. Path to HapMap3 SNP list (`w_hm3.snplist`).
#' @param chunksize Integer. Chunk size for munging (default 100000).
#' @param python Character. Python executable (default `"python"`).
#' @param tmp_dir Character. Optional output directory. If `NULL`, a new temp directory is created.
#' @param out_prefix Character. Optional output prefix within `tmp_dir` (default `"ldsc_gcim"`).
#'
#' @return A list with munging outputs.
#' @export
munge_ldsc_gcim <- function(munge_path,
                            glm_file,
                            interaction_term = NULL,
                            hm3_snplist,
                            chunksize = 100000,
                            python = "python",
                            tmp_dir = NULL,
                            out_prefix = "ldsc_gcim") {

  # ---- Helper functions ----
  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # ---- Resolve inputs if lists are provided ----
  gweis_source <- NULL
  
  if (is.list(glm_file)) {
    # Extract metadata
    if (!is.null(glm_file$meta)) {
      gweis_source <- glm_file$meta$trait_type
      # Auto-extract interaction_term if not provided
      if (is.null(interaction_term) && !is.null(glm_file$meta$int_term)) {
        interaction_term <- glm_file$meta$int_term
        message("Auto-detected interaction term: ", interaction_term)
      }
    }
    
    # Extract file path
    if (!is.null(glm_file$glm_file)) {
      glm_file <- glm_file$glm_file
    } else if (!is.null(glm_file$files$glm)) {
      glm_file <- glm_file$files$glm
    } else {
      .stop2("Input list does not contain 'glm_file' or 'files$glm' element.")
    }
  }

  # ---- Input validation ----
  if (!is.character(munge_path) || length(munge_path) != 1 || !nzchar(munge_path)) {
    .stop2("`munge_path` must be a non-empty character scalar.")
  }
  if (!file.exists(munge_path)) {
    .stop2("munge_sumstats.py not found at: ", munge_path)
  }

  if (!is.character(glm_file) || length(glm_file) != 1 || !nzchar(glm_file)) {
    .stop2("`glm_file` must be a non-empty character scalar path.")
  }
  if (!file.exists(glm_file)) {
    .stop2("PLINK GLM file not found: ", glm_file)
  }

  if (is.null(interaction_term) || !is.character(interaction_term) || 
      length(interaction_term) != 1 || !nzchar(interaction_term)) {
    .stop2(
      "`interaction_term` must be a non-empty string (e.g., 'ADDxCOVAR1').",
      "\nIf passing a list from q_gweis/b_gweis, the term should be auto-extracted."
    )
  }

  if (!is.character(hm3_snplist) || length(hm3_snplist) != 1 || !nzchar(hm3_snplist)) {
    .stop2("`hm3_snplist` must be a non-empty character scalar path.")
  }
  if (!file.exists(hm3_snplist)) {
    .stop2("HapMap3 SNP list not found: ", hm3_snplist)
  }

  if (!is.numeric(chunksize) || length(chunksize) != 1 || is.na(chunksize) || chunksize < 1000) {
    .stop2("`chunksize` must be a positive integer >= 1000")
  }
  chunksize <- as.integer(chunksize)

  if (!is.character(python) || length(python) != 1 || !nzchar(python)) {
    .stop2("`python` must be a non-empty character scalar.")
  }

  # ---- Create temporary directory ----
  if (is.null(tmp_dir)) {
    tmp_dir <- file.path(
      tempdir(),
      paste0("ldsc_gcim_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  }
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_prefix_path <- file.path(tmp_dir, out_prefix)
  munge_log <- paste0(out_prefix_path, "_munge.log")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 5 - Munge Summary Statistics for LDSC")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Interaction term: ", interaction_term)
  if (!is.null(gweis_source)) {
    message("GWEIS source: ", gweis_source, " trait")
  }

  # ---- Read GLM file ----
  glm <- tryCatch(
    utils::read.table(
      glm_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      comment.char = "",
      check.names = FALSE
    ),
    error = function(e) .stop2("Failed to read GLM file: ", glm_file, "\n", e$message)
  )

  if (nrow(glm) == 0) {
    .stop2("GLM file has 0 rows: ", glm_file)
  }

  if (!"TEST" %in% names(glm)) {
    .stop2("Column 'TEST' not found in GLM file: ", glm_file)
  }

  # ---- Extract interaction term ----
  sub <- glm[glm$TEST == interaction_term, , drop = FALSE]
  
  if (nrow(sub) == 0) {
    available_tests <- unique(glm$TEST)
    interaction_tests <- grep("^ADDxCOVAR", available_tests, value = TRUE)
    .stop2(
      "No rows found for interaction_term = '", interaction_term, "'.\n",
      "Available TEST values: ", paste(utils::head(available_tests, 20), collapse = ", "), "\n",
      "Available interaction terms: ",
      if (length(interaction_tests) > 0) paste(interaction_tests, collapse = ", ") else "none"
    )
  }

  message("SNPs with interaction term: ", nrow(sub))

  # ---- Check required columns ----
  core_needed <- c("ID", "A1", "P", "OBS_CT")
  miss_core <- setdiff(core_needed, names(sub))
  if (length(miss_core) > 0) {
    .stop2(
      "Missing required columns in GLM subset: ", paste(miss_core, collapse = ", "),
      "\nColumns available: ", paste(names(sub), collapse = ", ")
    )
  }

  # ---- Identify A2 column ----
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(sub)) {
      a2_col <- cand
      break
    }
  }
  if (is.null(a2_col)) {
    .stop2(
      "Could not find A2 column (tried A2/REF/ALT/ALT1).",
      "\nColumns available: ", paste(names(sub), collapse = ", ")
    )
  }

  message("Using allele 2 column: ", a2_col)

  # ---- Compute Z-scores ----
  z_source <- NULL
  
  if ("Z_STAT" %in% names(sub)) {
    z_source <- "Z_STAT"
    z <- suppressWarnings(as.numeric(sub$Z_STAT))
  } else if ("T_STAT" %in% names(sub)) {
    z_source <- "T_STAT"
    z <- suppressWarnings(as.numeric(sub$T_STAT))
  } else if (all(c("BETA", "SE") %in% names(sub))) {
    z_source <- "BETA/SE"
    beta <- suppressWarnings(as.numeric(sub$BETA))
    se   <- suppressWarnings(as.numeric(sub$SE))
    z <- beta / se
  } else {
    .stop2(
      "Could not find Z_STAT or T_STAT, and cannot compute Z (need BETA and SE).",
      "\nColumns available: ", paste(names(sub), collapse = ", ")
    )
  }

  message("Using Z-score from: ", z_source)

  # ---- Extract sample size ----
  n <- suppressWarnings(as.numeric(sub$OBS_CT))

  # ---- Build LDSC-ready data frame ----
  ldsc_df <- data.frame(
    SNP = as.character(sub$ID),
    A1  = as.character(sub$A1),
    A2  = as.character(sub[[a2_col]]),
    N   = n,
    Z   = z,
    stringsAsFactors = FALSE
  )

  # ---- Filter invalid entries ----
  keep <- stats::complete.cases(ldsc_df) & 
          is.finite(ldsc_df$Z) & 
          is.finite(ldsc_df$N) & 
          ldsc_df$N > 0
  
  n_removed <- sum(!keep)
  if (!all(keep)) {
    warning("Dropping ", n_removed, " SNPs with missing/invalid N or Z.", call. = FALSE)
    ldsc_df <- ldsc_df[keep, , drop = FALSE]
  }
  
  if (nrow(ldsc_df) == 0) {
    .stop2("No valid SNPs remain after filtering invalid N/Z.")
  }

  # ---- Remove duplicate SNP IDs ----
  n_duplicates <- sum(duplicated(ldsc_df$SNP))
  if (n_duplicates > 0) {
    warning(
      "Found ", n_duplicates, " duplicated SNP IDs; keeping first occurrence for each SNP.",
      call. = FALSE
    )
    ldsc_df <- ldsc_df[!duplicated(ldsc_df$SNP), , drop = FALSE]
  }

  n_snps <- nrow(ldsc_df)
  message("Valid SNPs for munging: ", n_snps)

  # ---- Write LDSC-ready summary statistics ----
  ldsc_sumstats <- file.path(tmp_dir, "ldsc_ready_sumstats.tsv")
  utils::write.table(
    ldsc_df,
    file = ldsc_sumstats,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  if (!file.exists(ldsc_sumstats)) {
    .stop2("Failed to create LDSC-ready sumstats file: ", ldsc_sumstats)
  }

  # ---- Run munge_sumstats.py ----
  message("Running munge_sumstats.py...")
  
  munge_args <- c(
    munge_path,
    "--sumstats", ldsc_sumstats,
    "--snp", "SNP",
    "--a1", "A1",
    "--a2", "A2",
    "--N-col", "N",
    "--signed-sumstats", "Z,0",
    "--chunksize", as.character(chunksize),
    "--merge-alleles", hm3_snplist,
    "--out", out_prefix_path
  )

  out <- tryCatch(
    system2(python, args = munge_args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      writeLines(paste("ERROR:", e$message), con = munge_log)
      structure(character(0), status = 1)
    }
  )
  
  writeLines(out, con = munge_log)
  
  status <- attr(out, "status")
  if (is.null(status)) status <- 0
  
  if (status != 0) {
    .stop2(
      "munge_sumstats.py exited with non-zero status (", status, ").",
      "\nCheck munge log at: ", munge_log
    )
  }

  # ---- Verify munged output ----
  munged_sumstats <- paste0(out_prefix_path, ".sumstats.gz")
  if (!file.exists(munged_sumstats)) {
    .stop2(
      "Munged sumstats file not found: ", munged_sumstats,
      "\nCheck munge log at: ", munge_log
    )
  }

  message(.bar())
  message("Munging completed successfully!")
  message("  Input SNPs: ", nrow(sub))
  message("  Valid SNPs: ", n_snps)
  if (n_removed > 0) {
    message("  SNPs removed (invalid): ", n_removed)
  }
  if (n_duplicates > 0) {
    message("  SNPs removed (duplicates): ", n_duplicates)
  }
  message("  LDSC-ready file: ", ldsc_sumstats)
  message("  Munged file: ", munged_sumstats)
  message("  Munge log: ", munge_log)
  message(.bar())

  # ---- Return results (pipeline-friendly structure) ----
  res <- list(
    # Backwards compatible fields
    ldsc_sumstats = ldsc_sumstats,
    munged_sumstats = munged_sumstats,
    munge_log = munge_log,
    tmp_dir = tmp_dir,
    n_snps = n_snps,

    # Pipeline-friendly additions
    files = list(
      ldsc_ready = ldsc_sumstats,
      munged = munged_sumstats,
      munge_log = munge_log
    ),
    meta = list(
      step = "MUNGE",
      interaction_term = interaction_term,
      gweis_source = gweis_source,
      a2_col = a2_col,
      z_source = z_source,
      chunksize = chunksize,
      n_snps_input = nrow(sub),
      n_snps_valid = n_snps,
      n_snps_removed = n_removed,
      n_snps_duplicates = n_duplicates,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}