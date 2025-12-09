#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(MetBrewer)
  library(scales)
  library(patchwork)
})

# ------------------------------ CLI ------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript plot_run_summary.R --final <final_qc_with_ashleys.tsv> --out <run_summary.pdf>\n",
      "Optional thresholds as env vars (override defaults):\n",
      "  TH_MAP_LT=0.5  TH_DUP_GT=0.7  TH_ERR_GT=0.01  TH_CV_GT=1.0  TH_GINI_GT=0.4\n",
      "  TH_GCABS_GT=0.6 TH_FRACCOV_LT=0.3 TH_DEPTH_LT=5e5 TH_SAT_LT=0.6 TH_ASH_P_GT=0.8\n", sep = "")
  quit(status = 1)
}
get_arg <- function(key, default = NULL) {
  pos <- which(args == key)
  if (length(pos) == 1 && pos < length(args)) args[pos + 1] else default
}
final_path <- get_arg("--final")
out_pdf    <- get_arg("--out", "plots/run_summary.pdf")

if (is.null(final_path)) {
  stop("Missing --final <path>")
}

# --------------------------- Params & Palettes -------------------------------
# Thresholds (overridable via env)
TH_MAP_LT     <- as.numeric(Sys.getenv("TH_MAP_LT",     "0.5"))
TH_DUP_GT     <- as.numeric(Sys.getenv("TH_DUP_GT",     "0.7"))
TH_ERR_GT     <- as.numeric(Sys.getenv("TH_ERR_GT",     "0.01"))
TH_CV_GT      <- as.numeric(Sys.getenv("TH_CV_GT",      "1.0"))
TH_GINI_GT    <- as.numeric(Sys.getenv("TH_GINI_GT",    "0.4"))
TH_GCABS_GT   <- as.numeric(Sys.getenv("TH_GCABS_GT",   "0.6"))
TH_FRACCOV_LT <- as.numeric(Sys.getenv("TH_FRACCOV_LT", "0.3"))
TH_DEPTH_LT   <- as.numeric(Sys.getenv("TH_DEPTH_LT",   "500000"))
TH_SAT_LT     <- as.numeric(Sys.getenv("TH_SAT_LT",     "0.6"))
TH_ASH_P_GT   <- as.numeric(Sys.getenv("TH_ASH_P_GT",   "0.8"))

# Nice palettes
pal_disc <- function(n = 8, name = "Hiroshige") met.brewer(name, max(3, min(21, n)))
pal_cont <- function(name = "Hokusai3") scale_color_gradientn(colours = met.brewer(name, 9))
pal_fill_cont <- function(name = "Hokusai3") scale_fill_gradientn(colours = met.brewer(name, 9))

theme_set(theme_pubr(base_size = 10, border = TRUE))

# ----------------------------- Load Data -------------------------------------
dat <- read_tsv(final_path, show_col_types = FALSE, guess_max = 1e6)

# Defensive renames (some columns can come in with _x/_y)
rename_if_present <- function(df, from, to) {
  if (from %in% names(df)) df <- dplyr::rename(df, !!to := !!rlang::sym(from))
  df
}
dat <- dat %>%
  rename_if_present("MedianMAPQ_x", "MedianMAPQ") %>%
  rename_if_present("MedianInsertSize_x", "MedianInsertSize") %>%
  rename_if_present("MedianInsertSize_y", "MedianInsertSize_dup") %>%
  rename_if_present("MedianMAPQ_y", "MedianMAPQ_dup")

# Ensure key columns exist
needed <- c("Library","Sample","MappedFraction","UnmappedFraction","DuplicateFraction",
            "ErrorRate","MedianMAPQ","total.read.count","avg.read.count","FractionCovered",
            "coverage_cv","coverage_gini","fold80_penalty","gc_pearson_r",
            "preseq_saturation","prediction","probability")
for (k in needed) if (!k %in% names(dat)) dat[[k]] <- NA

# ----------------------------- QC Flags --------------------------------------
flag_rules <- list(
  LowMapping       = ~ MappedFraction < TH_MAP_LT,
  HighDuplicates   = ~ DuplicateFraction > TH_DUP_GT,
  HighError        = ~ ErrorRate > TH_ERR_GT,
  PoorUniformity   = ~ coverage_cv > TH_CV_GT | coverage_gini > TH_GINI_GT,
  StrongGCBias     = ~ abs(gc_pearson_r) > TH_GCABS_GT,
  LowCoverage      = ~ FractionCovered < TH_FRACCOV_LT | total.read.count < TH_DEPTH_LT,
  LowComplexity    = ~ preseq_saturation < TH_SAT_LT,
  AshleyRisk       = ~ (!is.na(prediction) & prediction == 1 & !is.na(probability) & probability > TH_ASH_P_GT)
)

apply_flags <- function(df) {
  .flags <- lapply(flag_rules, function(fn) {
    tryCatch(as.logical(fn), error = function(e) rep(FALSE, nrow(df)))
  })
  fmat <- as.data.frame(.flags)
  names(fmat) <- names(flag_rules)
  out <- bind_cols(df, fmat)
  out$QC_FlagCount <- rowSums(fmat, na.rm = TRUE)
  out$QC_Flags <- apply(fmat, 1, function(r) {
    labs <- names(which(r))
    if (length(labs) == 0) "" else paste(labs, collapse = ",")
  })
  out
}
dat <- apply_flags(dat)

# ----------------------------- KPIs ------------------------------------------
kpi <- list(
  n_libraries = nrow(dat),
  med_map     = median(dat$MappedFraction, na.rm = TRUE),
  med_dup     = median(dat$DuplicateFraction, na.rm = TRUE),
  med_err     = median(dat$ErrorRate, na.rm = TRUE),
  med_depth   = median(dat$total.read.count, na.rm = TRUE),
  flagged     = mean(dat$QC_FlagCount > 0, na.rm = TRUE)
)

kpi_df <- tibble::tibble(
  Metric = c("Libraries","Median mapped","Median duplicates","Median error",
             "Median reads","% flagged"),
  Value  = c(
    kpi$n_libraries,
    percent(kpi$med_map, 0.1),
    percent(kpi$med_dup, 0.1),
    percent(kpi$med_err, 0.01),
    comma(kpi$med_depth),
    percent(kpi$flagged, 0.1)
  )
)

kpi_plot <- ggplot(kpi_df, aes(Metric, 1, fill = Metric)) +
  geom_tile(color = "grey30", height = 0.9, width = 0.95) +
  geom_text(aes(label = Value), fontface = 2, size = 4) +
  scale_fill_manual(values = pal_disc(nrow(kpi_df))) +
  labs(title = "Run overview (KPIs)") +
  theme_pubr(border = TRUE) +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(face = "bold"))

# -------------------------- Distributions page -------------------------------
dist_specs <- list(
  Mapping    = c("MappedFraction","UnmappedFraction","DuplicateFraction","MedianMAPQ"),
  Coverage   = c("total.read.count","avg.read.count","FractionCovered","MedianCoverage"),
  Uniformity = c("coverage_cv","coverage_gini","fold80_penalty","entropy"),
  Bias       = c("gc_pearson_r"),
  Preseq     = c("preseq_saturation"),
  Ashley     = c("probability","prediction")
)

mk_density <- function(df, col, facet = TRUE) {
  d <- df %>% mutate(.x = .data[[col]])
  ggplot(d, aes(.x, colour = Sample, fill = Sample)) +
    geom_density(alpha = 0.15, na.rm = TRUE) +
    scale_color_manual(values = pal_disc(length(unique(d$Sample)))) +
    scale_fill_manual(values = pal_disc(length(unique(d$Sample)))) +
    labs(x = col, y = "Density", title = col) +
    theme_pubr(border = TRUE) +
    theme(legend.position = "none")
}

mk_hist <- function(df, col) {
  d <- df %>% mutate(.x = .data[[col]])
  ggplot(d, aes(.x, fill = Sample)) +
    geom_histogram(bins = 40, alpha = 0.9, position = "identity", colour = "grey25") +
    scale_fill_manual(values = pal_disc(length(unique(d$Sample)))) +
    labs(x = col, y = "Count", title = col) +
    theme_pubr(border = TRUE) +
    theme(legend.position = "none")
}

dist_plots <- list(
  wrap_plots(lapply(dist_specs$Mapping,  mk_density, df = dat), ncol = 2) + plot_annotation(title = "Mapping distributions"),
  wrap_plots(lapply(dist_specs$Coverage, mk_density, df = dat), ncol = 2) + plot_annotation(title = "Coverage distributions"),
  wrap_plots(lapply(dist_specs$Uniformity,function(x) if (x %in% names(dat)) mk_density(dat, x) else NULL), ncol = 2) + plot_annotation(title = "Uniformity distributions"),
  wrap_plots(lapply(dist_specs$Bias,     function(x) if (x %in% names(dat)) mk_density(dat, x) else NULL), ncol = 1) + plot_annotation(title = "GC bias"),
  wrap_plots(lapply(dist_specs$Preseq,   function(x) if (x %in% names(dat)) mk_density(dat, x) else NULL), ncol = 1) + plot_annotation(title = "Library complexity (preseq)"),
  wrap_plots(list(
    mk_density(dat, "probability"),
    mk_hist(dat, "prediction")
  ), ncol = 2) + plot_annotation(title = "Ashley predictions")
)

# --------------------------- Scatter panels ----------------------------------
scatter <- function(x, y, color = "QC_FlagCount", title = NULL) {
  ggplot(dat, aes(.data[[x]], .data[[y]], color = .data[[color]], shape = Sample)) +
    geom_point(alpha = 0.9, size = 2) +
    pal_cont() + scale_shape_manual(values = seq_len(length(unique(dat$Sample)))) +
    labs(x = x, y = y, title = title) +
    theme_pubr(border = TRUE) +
    theme(legend.position = "right")
}

scatters <- wrap_plots(
  list(
    scatter("total.read.count","DuplicateFraction","QC_FlagCount","Duplication vs Read Depth"),
    scatter("coverage_cv","gc_pearson_r","QC_FlagCount","Uniformity (CV) vs GC bias"),
    scatter("fold80_penalty","total.read.count","QC_FlagCount","Fold-80 vs Read Depth"),
    scatter("ErrorRate","MappedFraction","QC_FlagCount","Error vs Mapping"),
    scatter("preseq_saturation","DuplicateFraction","QC_FlagCount","Saturation vs Duplicates")
  ),
  ncol = 2
)

# --------------------------- Flag summary ------------------------------------
flag_cols <- names(flag_rules)
flag_long <- dat %>%
  select(Sample, all_of(flag_cols)) %>%
  pivot_longer(-Sample, names_to = "Flag", values_to = "Triggered") %>%
  group_by(Sample, Flag) %>%
  summarize(n = sum(Triggered, na.rm = TRUE), .groups = "drop")

flag_plot <- ggplot(flag_long, aes(Flag, n, fill = Sample)) +
  geom_col(position = position_stack(), colour = "grey20") +
  coord_flip() +
  scale_fill_manual(values = pal_disc(length(unique(flag_long$Sample)))) +
  labs(title = "Flag prevalence by Sample", x = "", y = "Libraries flagged") +
  theme_pubr(border = TRUE)

# --------------------------- Correlation heatmap -----------------------------
corr_cols <- c("MappedFraction","DuplicateFraction","ErrorRate","total.read.count",
               "coverage_cv","coverage_gini","fold80_penalty","gc_pearson_r",
               "preseq_saturation","probability")
corr_df <- dat %>% select(any_of(corr_cols)) %>% mutate(across(everything(), as.numeric))
C <- suppressWarnings(cor(corr_df, use = "pairwise.complete.obs", method = "spearman"))
Cdf <- as.data.frame(as.table(C))
names(Cdf) <- c("Var1","Var2","rho")

heat <- ggplot(Cdf, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  pal_fill_cont() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  labs(title = "Spearman correlation (selected metrics)", x = "", y = "") +
  theme_pubr(border = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --------------------------- Outlier table (CSV) -----------------------------
outliers <- dat %>%
  arrange(desc(QC_FlagCount)) %>%
  select(Library, Sample, QC_FlagCount, QC_Flags,
         MappedFraction, DuplicateFraction, ErrorRate, total.read.count,
         coverage_cv, coverage_gini, fold80_penalty, gc_pearson_r,
         preseq_saturation, prediction, probability) %>%
  head(20)

# Save a companion CSV for quick triage
out_dir <- dirname(out_pdf)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write_tsv(outliers, file.path(out_dir, "qc_top_outliers.tsv"))
write_tsv(dat %>% select(Library, Sample, QC_FlagCount, QC_Flags),
          file.path(out_dir, "qc_flags.tsv"))

# --------------------------- Build the PDF -----------------------------------
pdf(out_pdf, width = 11.7, height = 8.3)  # A4 landscape
# Page 1: KPIs + flag summary
print(
  (kpi_plot / flag_plot) +
    plot_layout(heights = c(1, 2)) +
    plot_annotation(title = "Library QC — Run Summary")
)

# Pages: distributions
for (p in dist_plots) print(p)

# Page: scatter panels
print(scatters + plot_annotation(title = "Key relationships"))

# Page: correlation heatmap
print(heat)

dev.off()

cat("[plot_run_summary] Wrote: ", out_pdf, "\n", sep = "")
cat("[plot_run_summary] Extras:\n  - ", file.path(out_dir, "qc_top_outliers.tsv"),
    "\n  - ", file.path(out_dir, "qc_flags.tsv"), "\n", sep = "")
