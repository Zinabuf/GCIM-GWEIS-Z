#' Munge interaction summary statistics for LDSC
#'
#' @importFrom utils head
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
#' @param interaction_term Character. Interaction test name, e.g. `"ADDxCOVAR1"`.
#' @param hm3_snplist Character. Path to HapMap3 SNP list (`w_hm3.snplist`).
#' @param chunksize Integer. Chunk size for munging.
#' @param python Character. Python executable (default `"python"`).
#' @param tmp_dir Character. Optional output directory. If `NULL`, a new temp directory is created.
#' @param out_prefix Character. Optional output prefix within `tmp_dir` (default `"ldsc_gcim"`).
#'
#' @return A list with munging outputs.
#' @export
munge_ldsc_gcim <- function(munge_path,
                            glm_file,
                            interaction_term,
                            hm3_snplist,
                            chunksize = 100000,
                            python = "python",
                            tmp_dir = NULL,
                            out_prefix = "ldsc_gcim") {

  .stop2 <- function(...) stop(paste0(...), call. = FALSE)

  # allow passing prior-step lists
  if (is.list(glm_file) && !is.null(glm_file$glm_file)) {
    glm_file <- glm_file$glm_file
  }
  if (is.list(glm_file) && !is.null(glm_file$files$glm)) {
    glm_file <- glm_file$files$glm
  }

  if (!file.exists(munge_path)) .stop2("munge_sumstats.py not found at: ", munge_path)
  if (!is.character(glm_file) || length(glm_file) != 1 || !file.exists(glm_file)) .stop2("PLINK GLM file not found: ", glm_file)
  if (!file.exists(hm3_snplist)) .stop2("HapMap3 SNP list not found: ", hm3_snplist)

  if (!is.character(interaction_term) || length(interaction_term) != 1 || !nzchar(interaction_term)) {
    .stop2("interaction_term must be a non-empty string (e.g., 'ADDxCOVAR1').")
  }

  if (!is.numeric(chunksize) || length(chunksize) != 1 || is.na(chunksize) || chunksize < 1000) {
    .stop2("chunksize must be a positive integer >= 1000")
  }
  chunksize <- as.integer(chunksize)

  if (is.null(tmp_dir)) {
    tmp_dir <- file.path(tempdir(), paste0("ldsc_gcim_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix_path <- file.path(tmp_dir, out_prefix)

  munge_log <- paste0(out_prefix_path, "_munge.log")

  # ---- Read GLM ----
  glm <- tryCatch(
    utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE),
    error = function(e) .stop2("Failed to read GLM file: ", glm_file, "\n", e$message)
  )
  if (!"TEST" %in% names(glm)) .stop2("Column 'TEST' not found in GLM file: ", glm_file)

  sub <- glm[glm$TEST == interaction_term, , drop = FALSE]
  if (nrow(sub) == 0) {
    .stop2(
      "No rows found for interaction_term = ", interaction_term, "\n",
      "Example TEST values: ", paste(head(unique(glm$TEST), 20), collapse = ", ")
    )
  }

  core_needed <- c("ID", "A1", "P", "OBS_CT")
  miss_core <- setdiff(core_needed, names(sub))
  if (length(miss_core) > 0) {
    .stop2(
      "Missing required columns in GLM subset: ", paste(miss_core, collapse = ", "),
      "\nColumns available: ", paste(names(sub), collapse = ", ")
    )
  }

  # A2 column
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(sub)) { a2_col <- cand; break }
  }
  if (is.null(a2_col)) {
    .stop2("Could not find A2 column (tried A2/REF/ALT/ALT1). Columns: ", paste(names(sub), collapse = ", "))
  }

  # Z column: prefer Z_STAT; else T_STAT; else BETA/SE
  if ("Z_STAT" %in% names(sub)) {
    z <- suppressWarnings(as.numeric(sub$Z_STAT))
  } else if ("T_STAT" %in% names(sub)) {
    z <- suppressWarnings(as.numeric(sub$T_STAT))
  } else if (all(c("BETA", "SE") %in% names(sub))) {
    beta <- suppressWarnings(as.numeric(sub$BETA))
    se   <- suppressWarnings(as.numeric(sub$SE))
    z <- beta / se
  } else {
    .stop2("Could not find Z_STAT or T_STAT, and cannot compute Z (need BETA and SE). Columns: ",
           paste(names(sub), collapse = ", "))
  }

  n <- suppressWarnings(as.numeric(sub$OBS_CT))

  ldsc_df <- data.frame(
    SNP = as.character(sub$ID),
    A1  = as.character(sub$A1),
    A2  = as.character(sub[[a2_col]]),
    N   = n,
    Z   = z,
    stringsAsFactors = FALSE
  )

  keep <- stats::complete.cases(ldsc_df) & is.finite(ldsc_df$Z) & is.finite(ldsc_df$N) & ldsc_df$N > 0
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " SNPs with missing/invalid N or Z.", call. = FALSE)
    ldsc_df <- ldsc_df[keep, , drop = FALSE]
  }
  if (nrow(ldsc_df) == 0) .stop2("No valid SNPs remain after filtering invalid N/Z.")

  # Drop duplicate SNP IDs (munge_sumstats can fail otherwise)
  if (anyDuplicated(ldsc_df$SNP)) {
    warning("Found duplicated SNP IDs; keeping first occurrence for each SNP.", call. = FALSE)
    ldsc_df <- ldsc_df[!duplicated(ldsc_df$SNP), , drop = FALSE]
  }

  ldsc_sumstats <- file.path(tmp_dir, "ldsc_ready_sumstats.tsv")
  utils::write.table(
    ldsc_df,
    file = ldsc_sumstats,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )
  n_snps <- nrow(ldsc_df)

  # ---- Run munge_sumstats.py ----
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
  status <- attr(out, "status"); if (is.null(status)) status <- 0
  if (status != 0) {
    .stop2("munge_sumstats.py exited with non-zero status (", status, "). Check log: ", munge_log)
  }

  munged_sumstats <- paste0(out_prefix_path, ".sumstats.gz")
  if (!file.exists(munged_sumstats)) {
    .stop2("Munged sumstats file not found: ", munged_sumstats, "\nCheck munge log at: ", munge_log)
  }

  invisible(list(
    ldsc_sumstats = ldsc_sumstats,
    munged_sumstats = munged_sumstats,
    munge_log = munge_log,
    tmp_dir = tmp_dir,
    n_snps = n_snps,
    files = list(
      ldsc_ready = ldsc_sumstats,
      munged = munged_sumstats,
      munge_log = munge_log
    ),
    meta = list(
      step = "MUNGE",
      interaction_term = interaction_term,
      a2_col = a2_col,
      z_source = if ("Z_STAT" %in% names(sub)) "Z_STAT" else if ("T_STAT" %in% names(sub)) "T_STAT" else "BETA/SE",
      chunksize = chunksize
    )
  ))
}