#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

# ------------------------------- Helpers ------------------------------------

strip_bam <- function(x) {
  x <- basename(x)
  x <- str_replace(x, "\\.sort\\.mdup\\.bam$", "")
  x <- str_replace(x, "\\.sort\\.mdup$", "")
  x <- str_replace(x, "\\.mdup\\.bam$", "")
  x <- str_replace(x, "\\.mdup$", "")
  x <- str_replace(x, "\\.sort\\.bam$", "")
  x <- str_replace(x, "\\.sort$", "")
  x <- str_replace(x, "\\.bam$", "")
  x
}

# Find a column in df by case-insensitive match
find_col_ci <- function(df, candidates) {
  nm <- names(df)
  low <- tolower(nm)
  for (cand in candidates) {
    hit <- which(low == tolower(cand))
    if (length(hit) > 0) return(nm[hit[1]])
  }
  return(NULL)
}

# Remove a literal prefix if present (vectorized)
remove_literal_prefix <- function(x, prefix) {
  ok <- !is.na(prefix) & startsWith(x, prefix)
  x2 <- x
  x2[ok] <- substr(x[ok], nchar(prefix[ok]) + 1, nchar(x[ok]))
  x2
}


# robust reader
read_tsv_quiet <- function(path, ...) {
  readr::read_tsv(path, guess_max = 1e6, show_col_types = FALSE, progress = FALSE, ...)
}

# Build Library key in prediction.tsv: Library = paste0(sample, "_", strip_bam(cell))
normalize_pred_keys <- function(pred) {
  need <- c("cell","sample")
  miss <- setdiff(need, names(pred))
  if (length(miss) > 0) {
    stop("prediction.tsv missing required columns: ", paste(miss, collapse=", "))
  }

  pred %>%
    mutate(
      Base    = strip_bam(.data$cell),                 # A5573_L1_i301
      Library = paste0(.data$sample, "_", .data$Base, ".sort.mdup")
    )
}

# Add Base key to features: prefer 'sample_name', else 'cell'
normalize_feat_keys <- function(feat) {
  nm <- names(feat)
  key <- dplyr::case_when(
    "sample_name" %in% nm ~ "sample_name",
    "cell"        %in% nm ~ "cell",
    TRUE ~ NA_character_
  )
  if (is.na(key)) {
    warning("[merge_ashleys] features.tsv lacks 'sample_name' or 'cell', will skip merging features.")
    return(NULL)
  }
  feat %>% mutate(Base = strip_bam(.data[[key]]))
}

# Read a preseq curve for a Library (if present) and coerce key cols numeric
read_preseq_curve <- function(preseq_dir, lib) {
  f <- file.path(preseq_dir, paste0(lib, ".lc.tsv"))
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(read_tsv_quiet(f), error = function(e) NULL)
  if (is.null(df)) return(NULL)

  # Be tolerant to column types; coerce and keep only needed cols if present
  need <- c("TOTAL_READS","EXPECTED_DISTINCT")
  if (!all(need %in% names(df))) return(NULL)

  df <- df %>%
    mutate(
      TOTAL_READS       = suppressWarnings(as.numeric(.data$TOTAL_READS)),
      EXPECTED_DISTINCT = suppressWarnings(as.numeric(.data$EXPECTED_DISTINCT))
    )

  # Drop completely NA rows on both required columns
  df <- df %>% filter(!(is.na(TOTAL_READS) & is.na(EXPECTED_DISTINCT)))
  if (nrow(df) == 0) return(NULL)
  df
}

# Given a curve df and observed total reads, pick nearest point (no tidyr needed)
preseq_metrics_from_curve <- function(cur, observed_reads) {
  na_metrics <- tibble(
    preseq_distinct_at_observed = NA_real_,
    preseq_saturation           = NA_real_
  )

  if (is.null(cur) || !is.finite(observed_reads) || observed_reads <= 0) {
    return(na_metrics)
  }

  # If both columns are zero/NA everywhere, treat as degenerate
  all_zero_or_na <- {
    tr  <- cur$TOTAL_READS
    ed  <- cur$EXPECTED_DISTINCT
    # zeros or NAs only?
    (all(is.na(tr) | tr == 0) && all(is.na(ed) | ed == 0))
  }
  if (all_zero_or_na) return(na_metrics)

  # Find nearest row by TOTAL_READS (ignore NA)
  cur_ok <- cur %>% filter(is.finite(TOTAL_READS))
  if (nrow(cur_ok) == 0) return(na_metrics)

  idx <- which.min(abs(cur_ok$TOTAL_READS - observed_reads))
  distinct <- suppressWarnings(as.numeric(cur_ok$EXPECTED_DISTINCT[idx]))

  if (!is.finite(distinct) || observed_reads <= 0) return(na_metrics)

  tibble(
    preseq_distinct_at_observed = distinct,
    preseq_saturation           = distinct / observed_reads
  )
}


# Vectorized wrapper over a final table
attach_preseq_metrics <- function(final_df, out_root) {
  # Expect a numeric total.read.count and a Library ID
  if (!("Library" %in% names(final_df))) {
    warning("[merge_ashleys] final table lacks 'Library'; skipping preseq enrichment.")
    return(final_df)
  }
  trc <- "total.read.count"
  if (!(trc %in% names(final_df))) {
    warning("[merge_ashleys] final table lacks 'total.read.count'; skipping preseq enrichment.")
    return(final_df)
  }
  preseq_dir <- file.path(out_root, "preseq")

  enrich_row <- function(lib, obs) {
    cur <- read_preseq_curve(preseq_dir, lib)
    preseq_metrics_from_curve(cur, obs)
  }

  mets <- map2_dfr(final_df$Library, suppressWarnings(as.numeric(final_df[[trc]])), enrich_row)
  bind_cols(final_df, mets)
}

# ----------------------------- CLI & I/O -------------------------------------

ap <- argparse::ArgumentParser(description = "Merge final QC with Ashley prediction (and optional features), and enrich with preseq metrics.")
ap$add_argument("--final", required = TRUE, help = "Path to final_qc.tsv")
ap$add_argument("--pred",  required = TRUE, help = "Path to ashleys/prediction/prediction.tsv")
ap$add_argument("--feat",  required = FALSE, default = NULL, help = "Optional path to ashleys/features.tsv")
ap$add_argument("--out",   required = TRUE, help = "Output: final_qc_with_ashleys.tsv")
args <- ap$parse_args()

final_path <- normalizePath(args$final, mustWork = TRUE)
pred_path  <- normalizePath(args$pred,  mustWork = TRUE)
feat_path  <- if (!is.null(args$feat)) normalizePath(args$feat, mustWork = FALSE) else NULL
out_path   <- args$out

# OUT_ROOT = dirname(final_qc.tsv)
OUT_ROOT <- dirname(final_path)

message("[merge_ashleys] Reading inputs ...")
final <- read_tsv_quiet(final_path)
pred  <- read_tsv_quiet(pred_path)

# ---- Normalize final keys (Base, Sample, and a joinable Library key) ----
lib_col <- find_col_ci(final, c("Library"))
if (is.null(lib_col)) stop("[merge_ashleys] final_qc.tsv missing 'Library' column")

samp_col <- find_col_ci(final, c("sample", "Sample", "SAMPLE"))

final <- final %>%
  mutate(
    Library_raw = .data[[lib_col]],
    Sample_raw  = if (!is.null(samp_col)) .data[[samp_col]] else NA_character_
  )

# Compute Base
if (!is.null(samp_col)) {
  prefix <- paste0(final$Sample_raw, "_")
  lib_no_prefix <- remove_literal_prefix(final$Library_raw, prefix)
  final <- final %>% mutate(Base = strip_bam(lib_no_prefix))
} else {
  final <- final %>% mutate(Base = str_extract(.data$Library_raw, "A\\d+_L\\d+_i\\d+"))
}

# Compute Sample if missing
if (is.null(samp_col)) {
  final <- final %>% mutate(
    Sample_raw = str_replace(.data$Library_raw, "_A\\d+_L\\d+_i\\d+.*$", "")
  )
}

# Build normalized join key
final <- final %>% mutate(
  Library_key = paste0(.data$Sample_raw, "_", .data$Base, ".sort.mdup")
)

# ---- Normalize prediction keys to match final ----
pred <- normalize_pred_keys(pred) %>%
  mutate(Library_key = .data$Library)

message("[merge_ashleys] Key examples:")
message("  final Library_raw: ", final$Library_raw[1])
message("  final Sample_raw : ", final$Sample_raw[1])
message("  final Base       : ", final$Base[1])
message("  final Library_key: ", final$Library_key[1])
message("  pred  Library_key: ", pred$Library_key[1])

# ---- Join predictions (ONLY ONCE) ----
message("[merge_ashleys] Joining Ashley predictions ...")
keep_pred <- c("Library_key","prediction","probability","sample","cell")
pred_min  <- pred %>% select(any_of(keep_pred))
merged    <- final %>% left_join(pred_min, by = "Library_key")

message("[merge_ashleys] Matched predictions: ", sum(!is.na(merged$prediction)), " / ", nrow(merged))

# ---- Join Ashley features (by Base) ----
if (!is.null(feat_path) && file.exists(feat_path) && file.info(feat_path)$size > 0) {
  message("[merge_ashleys] Joining Ashley features ...")
  feat <- read_tsv_quiet(feat_path)

  feat <- normalize_feat_keys(feat)  # adds Base from sample_name/cell via strip_bam()
  if (!is.null(feat)) {
    # Prefix all feature columns except the join key to avoid name collisions
    feat_pref <- feat %>%
      # keep Base unprefixed
      rename_with(.cols = setdiff(names(feat), c("Base")), .fn = ~ paste0("feat_", .x))

    merged <- merged %>% left_join(feat_pref, by = "Base")

    # quick sanity: count how many rows got any non-NA feature value
    feat_cols <- grep("^feat_", names(merged), value = TRUE)
    if (length(feat_cols) > 0) {
      n_feat_rows <- sum(rowSums(!is.na(merged[, feat_cols, drop = FALSE])) > 0)
      message("[merge_ashleys] Rows with >=1 feature value: ", n_feat_rows, " / ", nrow(merged))
    }
  } else {
    message("[merge_ashleys] Skipping features merge (no usable key).")
  }
} else {
  message("[merge_ashleys] No features file supplied or empty; skipping features merge.")
}

# [preseq skipped] Do not enrich with preseq metrics for now.
message("[merge_ashleys] Skipping preseq enrichment (requested).")

# Robust/atomic write: write to tmp, verify, then rename to final
out_path <- normalizePath(args$out, mustWork = FALSE)
out_dir  <- dirname(out_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# add a little context into the log
message("[merge_ashleys] Writing (atomic) to: ", out_path)
message("[merge_ashleys] getwd(): ", getwd())

tmp_path <- paste0(out_path, ".tmp")

# write tmp
readr::write_tsv(merged, tmp_path)

# verify tmp exists and is non-empty
if (!file.exists(tmp_path)) {
  stop("[merge_ashleys] ERROR: tmp not found after write: ", tmp_path)
}
sz <- tryCatch(file.size(tmp_path), error = function(e) NA_real_)
if (!is.finite(sz) || sz <= 0) {
  stop("[merge_ashleys] ERROR: tmp file is empty: ", tmp_path)
}

# move tmp -> final (atomic on same filesystem)
ok <- file.rename(tmp_path, out_path)
if (!isTRUE(ok)) {
  # fallback for weird FS semantics
  file.copy(tmp_path, out_path, overwrite = TRUE)
  unlink(tmp_path)
}

# verify final exists and is non-empty
if (!file.exists(out_path)) {
  stop("[merge_ashleys] ERROR: final not found after write: ", out_path)
}
final_sz <- tryCatch(file.size(out_path), error = function(e) NA_real_)
if (!is.finite(final_sz) || final_sz <= 0) {
  stop("[merge_ashleys] ERROR: final file is empty: ", out_path)
}

message("[merge_ashleys] Wrote OK (", final_sz, " bytes): ", out_path)
message("[merge_ashleys] Done.")