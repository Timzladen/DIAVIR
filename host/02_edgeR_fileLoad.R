#######################################
# edgeR analysis of MMTESP mapped reads
#######################################
#set path
base_path <- "/scratch/timtd/transcriptomes/MMETSP/fasta"

#load sample data
samples <- read.csv(
  "/scratch/timtd/transcriptomes/trimmed/salmon_quant/sample_data.csv",
  stringsAsFactors = TRUE
)

#load salmon quantifiactions
files <- file.path(base_path, "salmon_quant_renamed", as.character(samples$sample), "quant.sf")
names(files) <- as.character(samples$sample)

#load annotation table
eggnog_full_named <- read.delim(
  "/scratch/timtd/transcriptomes/MMETSP/fasta/deseq2_prereview/pn_eggnog.renamed.tsv",
  header = TRUE, sep = "\t", quote = ""
)

# ---- output ----
outdir <- file.path(base_path, "edger_anlysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


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


# Keep only IDs present in counts
eggnog_full_named <- eggnog_full_named[eggnog_full_named$X.query %in% rownames(counts), , drop = FALSE]

message("Annotation rows after cleaning + matching to counts: ", nrow(eggnog_full_named))


# ---- save outputs ----
saveRDS(samples, file.path(outdir, "sample_table.rds"))
saveRDS(counts, file.path(outdir, "counts_all.rds"))
saveRDS(txi, file.path(outdir, "tximport_all.rds"))
saveRDS(y_all, file.path(outdir, "edger_dge_all.rds"))
saveRDS(eggnog_full_named, file.path(outdir, "eggnog_full_named_matched.rds"))

