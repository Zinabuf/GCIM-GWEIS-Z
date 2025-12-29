#' Replace the third column with PRS
#'
#' Merges PRS values into a phenotype/covariate file and replaces column 3 with PRS.
#' Output is written to a temporary directory.
#'
#' @importFrom utils read.table write.table
#' @param dis_cov_file Character (path) OR list. Path to phenotype/covariate file
#'   (FID IID + >=1 columns). If a list is provided, it must contain `dis_cov_file`
#'   or `covar_file`.
#' @param prs_file Character (path) OR list from `prs_scores()` containing `$prs_file`.
#' @param out_name Character. Output filename (default: "tar_cov_prs.txt").
#' @param on_missing Character. What to do if PRS is missing for some samples:
#'   "stop" (default) or "drop".
#'
#' @return A list with:
#' \describe{
#'   \item{output_file}{Path to the output file with PRS replacing column 3}
#'   \item{tmp_dir}{Temporary directory used by the function}
#'   \item{n_samples}{Number of samples in the output}
#'   \item{n_cols}{Number of columns in the output}
#'   \item{col3_replaced}{Name of the original third column that was replaced}
#'   \item{files}{Named list of key output files (pipeline-friendly)}
#'   \item{meta}{Named list of metadata (pipeline-friendly)}
#' }
#'
#' @export
replace_covariate_with_prs <- function(dis_cov_file,
                                       prs_file,
                                       out_name = "tar_cov_prs.txt",
                                       on_missing = c("stop", "drop")) {

  # ---- Helper functions ----
  .bar <- function() paste(rep("=", 60), collapse = "")
  .stop2 <- function(...) stop(paste0(...), call. = FALSE)
  
  on_missing <- match.arg(on_missing)

  # ---- Allow passing previous-step lists ----
  original_cov_source <- NULL
  original_prs_source <- NULL
  
  if (is.list(dis_cov_file)) {
    if (!is.null(dis_cov_file$meta)) {
      original_cov_source <- dis_cov_file$meta$step
    }
    if (!is.null(dis_cov_file$dis_cov_file)) {
      dis_cov_file <- dis_cov_file$dis_cov_file
    } else if (!is.null(dis_cov_file$covar_file)) {
      dis_cov_file <- dis_cov_file$covar_file
    } else if (!is.null(dis_cov_file$output_file)) {
      dis_cov_file <- dis_cov_file$output_file
    } else {
      .stop2("Input list does not contain 'dis_cov_file', 'covar_file', or 'output_file'.")
    }
  }
  
  if (is.list(prs_file)) {
    if (!is.null(prs_file$meta)) {
      original_prs_source <- prs_file$meta$step
    }
    if (!is.null(prs_file$prs_file)) {
      prs_file <- prs_file$prs_file
    } else {
      .stop2("Input list does not contain 'prs_file' element.")
    }
  }

  # ---- Input validation ----
  if (!is.character(dis_cov_file) || length(dis_cov_file) != 1 || !nzchar(dis_cov_file)) {
    .stop2("`dis_cov_file` must be a non-empty character scalar path.")
  }
  if (!file.exists(dis_cov_file)) {
    .stop2("Covariate file not found: ", dis_cov_file)
  }

  if (!is.character(prs_file) || length(prs_file) != 1 || !nzchar(prs_file)) {
    .stop2("`prs_file` must be a non-empty character scalar path.")
  }
  if (!file.exists(prs_file)) {
    .stop2("PRS file not found: ", prs_file)
  }

  if (!is.character(out_name) || length(out_name) != 1 || !nzchar(out_name)) {
    .stop2("`out_name` must be a non-empty character scalar.")
  }

  # ---- Create function-specific temporary directory ----
  tmp_dir <- file.path(
    tempdir(),
    paste0("replace_cov_prs_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  message(.bar())
  message("GCIM-GWEIS-Z: Step 3 - Replace Exposure with PRS")
  message(.bar())
  message("Output directory: ", tmp_dir)

  # ---- Read covariate file ----
  covar <- tryCatch(
    utils::read.table(
      dis_cov_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read covariate file: ", dis_cov_file, "\n", e$message)
  )

  if (nrow(covar) == 0) {
    .stop2("Covariate file has 0 rows: ", dis_cov_file)
  }

  if (ncol(covar) < 3) {
    .stop2(
      "Covariate file must contain at least 3 columns (FID, IID, and column 3).",
      "\nActual columns: ", ncol(covar),
      "\nFile: ", dis_cov_file
    )
  }

  if (!all(c("FID", "IID") %in% names(covar))) {
    .stop2(
      "Covariate file must contain 'FID' and 'IID' columns.",
      "\nActual columns: ", paste(names(covar), collapse = ", "),
      "\nFile: ", dis_cov_file
    )
  }

  col3_name <- names(covar)[3]
  n_samples_original <- nrow(covar)
  
  message("Original covariate file:")
  message("  Samples: ", n_samples_original)
  message("  Columns: ", ncol(covar))
  message("  Column 3 to replace: '", col3_name, "'")

  # ---- Read PRS file ----
  prs <- tryCatch(
    utils::read.table(
      prs_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = ""
    ),
    error = function(e) .stop2("Failed to read PRS file: ", prs_file, "\n", e$message)
  )

  if (nrow(prs) == 0) {
    .stop2("PRS file has 0 rows: ", prs_file)
  }

  req <- c("FID", "IID", "PRS")
  miss <- setdiff(req, names(prs))
  if (length(miss) > 0) {
    .stop2(
      "PRS file missing required columns: ", paste(miss, collapse = ", "),
      "\nExpected columns: FID, IID, PRS",
      "\nActual columns: ", paste(names(prs), collapse = ", "),
      "\nFile: ", prs_file
    )
  }

  message("PRS file:")
  message("  Samples: ", nrow(prs))

  # ---- Check for duplicates in PRS ----
  key <- paste(prs$FID, prs$IID, sep = "___")
  if (anyDuplicated(key)) {
    dups <- prs[duplicated(key) | duplicated(key, fromLast = TRUE), c("FID", "IID")]
    dup_file <- file.path(tmp_dir, "prs_duplicates.txt")
    utils::write.table(
      dups,
      file = dup_file,
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t",
      quote = FALSE
    )
    .stop2(
      "PRS file has duplicated FID/IID pairs. Please deduplicate before merging.",
      "\nDuplicated samples saved to: ", dup_file
    )
  }

  # ---- Match PRS to covariate file (preserve row order) ----
  cov_key <- paste(covar$FID, covar$IID, sep = "___")
  prs_map <- prs$PRS
  names(prs_map) <- key
  prs_vec <- unname(prs_map[cov_key])

  n_missing <- sum(is.na(prs_vec))
  missing_file <- file.path(tmp_dir, "samples_missing_prs.txt")

  if (n_missing > 0) {
    missing_samples <- covar[is.na(prs_vec), c("FID", "IID"), drop = FALSE]
    utils::write.table(
      missing_samples,
      file = missing_file,
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t",
      quote = FALSE
    )

    if (on_missing == "stop") {
      .stop2(
        "Missing PRS values for ", n_missing, " samples (",
        sprintf("%.1f%%", 100 * n_missing / n_samples_original), ").\n",
        "Samples with missing PRS saved to: ", missing_file, "\n",
        "Suggestion: Check that FID/IID match between covariate and PRS files,\n",
        "  or rerun with on_missing='drop' to exclude these samples."
      )
    } else {
      warning(
        "Dropping ", n_missing, " samples (",
        sprintf("%.1f%%", 100 * n_missing / n_samples_original),
        ") with missing PRS.\n",
        "  See: ", missing_file,
        call. = FALSE
      )
      keep <- !is.na(prs_vec)
      covar <- covar[keep, , drop = FALSE]
      prs_vec <- prs_vec[keep]
    }
  }

  if (nrow(covar) == 0) {
    .stop2("No samples remaining after PRS matching.")
  }

  # ---- Replace column 3 with PRS ----
  covar[[3]] <- prs_vec
  names(covar)[3] <- "PRS"  # Rename column 3 to "PRS" for clarity

  # ---- Save output ----
  out_file <- file.path(tmp_dir, out_name)
  utils::write.table(
    covar,
    file = out_file,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
  )

  if (!file.exists(out_file)) {
    .stop2("Failed to create output file: ", out_file)
  }

  message(.bar())
  message("PRS replacement completed successfully!")
  message("  Original column 3: '", col3_name, "'")
  message("  New column 3: 'PRS'")
  message("  Output samples: ", nrow(covar), " / ", n_samples_original)
  if (n_missing > 0) {
    message("  Samples dropped: ", n_missing)
  }
  message("  Output file: ", out_file)
  message(.bar())

  # ---- Return results (pipeline-friendly structure) ----
  res <- list(
    # Backwards compatible fields
    output_file = out_file,
    tmp_dir = tmp_dir,
    n_samples = nrow(covar),
    n_cols = ncol(covar),
    col3_replaced = col3_name,

    # Pipeline-friendly additions
    files = list(
      covar_prs = out_file,
      missing_prs = if (file.exists(missing_file)) missing_file else NA_character_
    ),
    meta = list(
      step = "PRS_REPLACE",
      replaced_col = col3_name,
      new_col_name = "PRS",
      on_missing = on_missing,
      n_samples_original = n_samples_original,
      n_samples_final = nrow(covar),
      n_samples_dropped = n_missing,
      source_covar = original_cov_source,
      source_prs = original_prs_source,
      timestamp = Sys.time()
    )
  )

  invisible(res)
}