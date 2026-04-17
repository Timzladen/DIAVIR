#!/usr/bin/env Rscript
# =========================
# 01_build_edger_objects.R
# Build edgeR-ready objects from Salmon quant directories + eggNOG annotation
# =========================

# ---- paths you provided ----
base_path <- "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst"

samples <- read.csv(
  "/scratch/timtd/transcriptomes/trimmed/salmon_quant/sample_data.csv",
  stringsAsFactors = TRUE
)

files <- file.path(base_path, as.character(samples$sample), "quant.sf")
names(files) <- as.character(samples$sample)

eggnog_full_named <- read.delim(
  "/scratch/timtd/transcriptomes/trimmed/salmon_quant/eggnog_full_named_clean.tsv",
  header = TRUE, sep = "\t", quote = ""
)

# ---- output ----
outdir <- file.path(base_path, "edger_objects")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- install/load packages ----
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

install_if_missing(c("edgeR", "tximport", "stringr"))

suppressPackageStartupMessages({
  library(edgeR)
  library(tximport)
  library(stringr)
})

# ---- basic checks ----
if (!("sample" %in% colnames(samples))) {
  stop("samples metadata must contain a column named 'sample'")
}
if (!file.exists(files[1])) {
  missing_files <- files[!file.exists(files)]
  stop("Some quant.sf files do not exist. Example missing:\n  ",
       paste(head(missing_files, 10), collapse = "\n  "))
}

# ---- eggNOG: add genus from spname ----
if (!("spname" %in% colnames(eggnog_full_named))) {
  warning("Column 'spname' not found in eggnog_full_named; genus will not be added.")
} else {
  eggnog_full_named$genus <- word(eggnog_full_named$spname, 1)
}

# ---- import Salmon counts ----
txi <- tximport(
  files,
  type = "salmon",
  txOut = TRUE,
  countsFromAbundance = "no"
)

counts <- txi$counts

# ---- align samples metadata to counts columns ----
# Ensure same order and no missing
samples$sample <- as.character(samples$sample)
samples <- samples[match(colnames(counts), samples$sample), , drop = FALSE]

if (any(is.na(samples$sample))) {
  stop("Some samples in counts were not found in sample_data.csv. Missing:\n  ",
       paste(colnames(counts)[is.na(samples$sample)], collapse = ", "))
}

# ---- build edgeR object (all) ----
y_all <- DGEList(counts = counts)

# carry metadata columns into y$samples (except sample id itself)
meta_cols <- setdiff(colnames(samples), "sample")
if (length(meta_cols) > 0) {
  y_all$samples <- cbind(y_all$samples, samples[, meta_cols, drop = FALSE])
}

# ---- CLEAN eggNOG transcript_id (remove junk log lines) ----
if (!("transcript_id" %in% colnames(eggnog_full_named))) {
  stop("eggnog_full_named must contain column 'transcript_id'")
}

# Remove rows where transcript_id clearly isn't a contig/transcript ID
# Your IDs look like: NODE_..._g..._i...
is_real_id <- str_detect(eggnog_full_named$transcript_id, "^NODE_")
eggnog_full_named <- eggnog_full_named[is_real_id, , drop = FALSE]

# Keep only IDs present in counts
eggnog_full_named <- eggnog_full_named[eggnog_full_named$transcript_id %in% rownames(counts), , drop = FALSE]

message("Annotation rows after cleaning + matching to counts: ", nrow(eggnog_full_named))

# ---- define taxon_group using available columns ----
# We can do a robust bacteria call from domain.
# Diatom call needs clues, because domain alone won't distinguish diatoms from other eukaryotes.
eggnog_full_named$taxon_group <- NA_character_

if ("domain" %in% colnames(eggnog_full_named)) {
  dom <- tolower(as.character(eggnog_full_named$domain))
  eggnog_full_named$taxon_group[dom == "bacteria"] <- "bacteria"
}

# Diatom inference: only among eukaryotes (if domain exists), using organism_type and/or spname
# This is conservative; it won’t mislabel random eukaryotes as diatoms unless there are diatom keywords.
diatom_hits <- rep(FALSE, nrow(eggnog_full_named))

if ("organism_type" %in% colnames(eggnog_full_named)) {
  ot <- tolower(as.character(eggnog_full_named$organism_type))
  diatom_hits <- diatom_hits | str_detect(ot, "diatom|bacillariophy|pennate|centric")
}
if ("spname" %in% colnames(eggnog_full_named)) {
  sn <- tolower(as.character(eggnog_full_named$spname))
  # catches many eggNOG-style names when diatom is explicit
  diatom_hits <- diatom_hits | str_detect(sn, "bacillariophy|diatom")
}

# Apply diatom label; if domain column exists, restrict to eukaryota rows for safety
if ("domain" %in% colnames(eggnog_full_named)) {
  dom <- tolower(as.character(eggnog_full_named$domain))
  eggnog_full_named$taxon_group[dom == "eukaryota" & diatom_hits] <- "diatom"
} else {
  eggnog_full_named$taxon_group[diatom_hits] <- "diatom"
}

classified_n <- sum(!is.na(eggnog_full_named$taxon_group))
message("Taxon group assigned for ", classified_n, " / ", nrow(eggnog_full_named), " annotated contigs.")
message("  bacteria: ", sum(eggnog_full_named$taxon_group == "bacteria", na.rm = TRUE))
message("  diatom:   ", sum(eggnog_full_named$taxon_group == "diatom", na.rm = TRUE))

# ---- subset DGEList objects ----
subset_by_group <- function(y, annot_df, group_value) {
  keep_ids <- annot_df$transcript_id[annot_df$taxon_group == group_value]
  keep_ids <- intersect(keep_ids, rownames(y$counts))
  y[keep_ids, , keep.lib.sizes = FALSE]
}

y_bacteria <- NULL
y_diatom <- NULL

if (any(eggnog_full_named$taxon_group == "bacteria", na.rm = TRUE)) {
  y_bacteria <- subset_by_group(y_all, eggnog_full_named, "bacteria")
}
if (any(eggnog_full_named$taxon_group == "diatom", na.rm = TRUE)) {
  y_diatom <- subset_by_group(y_all, eggnog_full_named, "diatom")
}

# ---- save outputs ----
saveRDS(samples, file.path(outdir, "sample_table.rds"))
saveRDS(counts, file.path(outdir, "counts_all.rds"))
saveRDS(txi, file.path(outdir, "tximport_all.rds"))
saveRDS(y_all, file.path(outdir, "edger_dge_all.rds"))
saveRDS(eggnog_full_named, file.path(outdir, "eggnog_full_named_matched.rds"))

if (!is.null(y_bacteria) && nrow(y_bacteria$counts) > 0) {
  saveRDS(y_bacteria, file.path(outdir, "edger_dge_bacteria.rds"))
}
if (!is.null(y_diatom) && nrow(y_diatom$counts) > 0) {
  saveRDS(y_diatom, file.path(outdir, "edger_dge_diatom.rds"))
  saveRDS(y_diatom$counts, file.path(outdir, "counts_diatom.rds"))
}
