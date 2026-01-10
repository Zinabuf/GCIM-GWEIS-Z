#' Munge interaction summary statistics for LDSC
#'
#' Creates an LDSC-ready table (SNP, A1, A2, N, Z, P) from PLINK2 .glm output
#' and munges it with munge_sumstats.py.
#'
#' @importFrom utils read.table write.table tail
#' @importFrom stats complete.cases
#' @param munge_path Path to munge_sumstats.py.
#' @param glm_file Path to PLINK2 .glm file (or list from q_gweis/b_gweis).
#' @param interaction_term Interaction term, e.g. "ADDxCOVAR1".
#' @param hm3_snplist Path to HapMap3 SNP list.
#' @param chunksize Chunk size for munging large files.
#' @param python Python executable to run munge_sumstats.py.
#' @param tmp_dir Optional directory for outputs (default tempdir()).
#' @param out_prefix Output prefix for munged file(s).
#' @param show_munge_log_tail_on_error Logical; print log tail on failure.

#' @export
munge_ldsc_gcim <- function(munge_path,
                            glm_file,
                            interaction_term = NULL,
                            hm3_snplist,
                            chunksize = 100000,
                            python = "python",
                            tmp_dir = NULL,
                            out_prefix = "ldsc_gcim",
                            show_munge_log_tail_on_error = TRUE) {

  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)
  .tail_file <- function(path, n = 140) {
    if (!file.exists(path)) return(character(0))
    utils::tail(readLines(path, warn = FALSE), n)
  }

  # ---- Resolve list input (q_gweis/b_gweis output) ----
  gweis_source <- NULL
  if (is.list(glm_file)) {
    if (!is.null(glm_file$meta)) {
      gweis_source <- glm_file$meta$trait_type
      if (is.null(interaction_term) && !is.null(glm_file$meta$int_term)) {
        interaction_term <- glm_file$meta$int_term
        message("Auto-detected interaction term: ", interaction_term)
      }
    }
    if (!is.null(glm_file$glm_file)) glm_file <- glm_file$glm_file
    else if (!is.null(glm_file$files$glm)) glm_file <- glm_file$files$glm
    else .stop2("Input list does not contain 'glm_file' or 'files$glm'.")
  }

  # ---- Validate ----
  if (!file.exists(munge_path)) .stop2("munge_sumstats.py not found at: ", munge_path)
  if (!file.exists(glm_file)) .stop2("PLINK GLM file not found: ", glm_file)
  if (!file.exists(hm3_snplist)) .stop2("HapMap3 SNP list not found: ", hm3_snplist)

  chunksize <- as.integer(chunksize)
  if (is.na(chunksize) || chunksize < 1000) .stop2("`chunksize` must be an integer >= 1000")

  if (is.null(interaction_term) || !is.character(interaction_term) || !nzchar(interaction_term)) {
    .stop2("`interaction_term` must be provided (e.g., 'ADDxPRS') or auto-detected from q_gweis().")
  }

  # ---- Output dir ----
  if (is.null(tmp_dir)) {
    tmp_dir <- file.path(tempdir(), paste0("ldsc_gcim_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  out_prefix_path <- file.path(tmp_dir, out_prefix)
  munge_log <- paste0(out_prefix_path, "_munge.log")

  message(.bar())
  message("GCIM-GWEIS-Z: Step 5 - Munge Summary Statistics for LDSC")
  message(.bar())
  message("Output directory: ", tmp_dir)
  message("Interaction term: ", interaction_term)
  if (!is.null(gweis_source)) message("GWEIS source: ", gweis_source, " trait")

  # ---- Read GLM ----
  glm <- utils::read.table(glm_file, header = TRUE, stringsAsFactors = FALSE,
                           comment.char = "", check.names = FALSE)
  if (!"TEST" %in% names(glm)) .stop2("Column 'TEST' not found in: ", glm_file)

  sub <- glm[glm$TEST == interaction_term, , drop = FALSE]
  if (nrow(sub) == 0) .stop2("No rows found for interaction_term='", interaction_term, "' in: ", glm_file)

  message("SNPs with interaction term: ", nrow(sub))

  # ---- Required columns (NOTE: includes P) ----
  core_needed <- c("ID", "A1", "P", "OBS_CT")
  miss <- setdiff(core_needed, names(sub))
  if (length(miss) > 0) .stop2("Missing required columns in GLM subset: ", paste(miss, collapse = ", "))

  # ---- A2 column ----
  a2_col <- NULL
  for (cand in c("A2", "REF", "ALT", "ALT1")) {
    if (cand %in% names(sub)) { a2_col <- cand; break }
  }
  if (is.null(a2_col)) .stop2("Could not find allele2 column (A2/REF/ALT/ALT1).")
  message("Using allele 2 column: ", a2_col)

  # ---- Z ----
  if ("Z_STAT" %in% names(sub)) {
    z_source <- "Z_STAT"
    z <- suppressWarnings(as.numeric(sub$Z_STAT))
  } else if ("T_STAT" %in% names(sub)) {
    z_source <- "T_STAT"
    z <- suppressWarnings(as.numeric(sub$T_STAT))
  } else if (all(c("BETA", "SE") %in% names(sub))) {
    z_source <- "BETA/SE"
    z <- suppressWarnings(as.numeric(sub$BETA)) / suppressWarnings(as.numeric(sub$SE))
  } else {
    .stop2("Could not find Z_STAT or T_STAT, and cannot compute Z (need BETA and SE).")
  }
  message("Using Z-score from: ", z_source)

  n <- suppressWarnings(as.numeric(sub$OBS_CT))
  p <- suppressWarnings(as.numeric(sub$P))

  # ---- Build LDSC-ready df INCLUDING P (this fixes your error) ----
  ldsc_df <- data.frame(
    SNP = as.character(sub$ID),
    A1  = as.character(sub$A1),
    A2  = as.character(sub[[a2_col]]),
    N   = n,
    Z   = z,
    P   = p,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  keep <- stats::complete.cases(ldsc_df) &
    is.finite(ldsc_df$Z) &
    is.finite(ldsc_df$N) & ldsc_df$N > 0 &
    is.finite(ldsc_df$P) & ldsc_df$P > 0 & ldsc_df$P <= 1

  if (any(!keep)) {
    warning("Dropping ", sum(!keep), " SNPs with invalid/missing N/Z/P.", call. = FALSE)
    ldsc_df <- ldsc_df[keep, , drop = FALSE]
  }
  if (nrow(ldsc_df) == 0) .stop2("No valid SNPs remain after filtering invalid N/Z/P.")

  if (any(duplicated(ldsc_df$SNP))) {
    warning("Duplicated SNP IDs found; keeping first occurrence.", call. = FALSE)
    ldsc_df <- ldsc_df[!duplicated(ldsc_df$SNP), , drop = FALSE]
  }

  message("Valid SNPs for munging: ", nrow(ldsc_df))

  # ---- Write sumstats ----
  ldsc_sumstats <- file.path(tmp_dir, "ldsc_ready_sumstats.tsv")
  utils::write.table(ldsc_df, ldsc_sumstats, row.names = FALSE, col.names = TRUE,
                     sep = "\t", quote = FALSE)

  # ---- Run munge_sumstats.py (stdout+stderr -> log file) ----
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

  status <- system2(python, args = munge_args, stdout = munge_log, stderr = munge_log)
  if (is.null(status)) status <- 0

  if (status != 0) {
    if (isTRUE(show_munge_log_tail_on_error) && file.exists(munge_log)) {
      cat("\n---- munge_sumstats.py log (tail) ----\n")
      cat(paste(.tail_file(munge_log, n = 180), collapse = "\n"), "\n")
    }
    .stop2("munge_sumstats.py exited with non-zero status (", status, ").\nCheck munge log at: ", munge_log)
  }

  munged_sumstats <- paste0(out_prefix_path, ".sumstats.gz")
  if (!file.exists(munged_sumstats)) .stop2("Munged output not found: ", munged_sumstats)

  message(.bar())
  message("Munging completed successfully!")
  message("  LDSC-ready file: ", ldsc_sumstats)
  message("  Munged file: ", munged_sumstats)
  message("  Munge log: ", munge_log)
  message(.bar())

  invisible(list(
    ldsc_sumstats = ldsc_sumstats,
    munged_sumstats = munged_sumstats,
    munge_log = munge_log,
    tmp_dir = tmp_dir,
    meta = list(step = "MUNGE", interaction_term = interaction_term, z_source = z_source,
                a2_col = a2_col, chunksize = chunksize, timestamp = Sys.time())
  ))
}