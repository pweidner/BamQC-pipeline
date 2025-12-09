suppressPackageStartupMessages({
  library(optparse); library(readr); library(dplyr)
})

opt_list <- list(
  make_option(c("--final"), type="character"),
  make_option(c("--pred"),  type="character"),
  make_option(c("--feat"),  type="character", default=NA),
  make_option(c("--out"),   type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

final <- read_tsv(opt$final, show_col_types = FALSE)
pred  <- read_tsv(opt$pred,  show_col_types = FALSE)

# enforce Library presence
stopifnot("Library" %in% colnames(final))
stopifnot("Library" %in% colnames(pred))

# select small Ashley columns by default (prediction, probability, sample if present)
keep_pred <- intersect(c("Library","ashleys_prediction","ashleys_probability","ashleys_sample"), colnames(pred))
pred_slim <- pred %>% select(all_of(keep_pred)) %>% distinct()

merged <- final %>% left_join(pred_slim, by="Library")

# optional features
if (!is.na(opt$feat) && file.exists(opt$feat)) {
  feat <- read_tsv(opt$feat, show_col_types = FALSE)
  if (!"Library" %in% colnames(feat)) {
    warning("[merge_ashleys] features.norm.tsv lacks 'Library'; skipping.")
  } else {
    # join all feature columns; they’re already namespaced (e.g., W10_* etc.)
    merged <- merged %>% left_join(feat, by="Library")
  }
}

write_tsv(merged, opt$out)
