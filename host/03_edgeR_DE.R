#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

# =========================
# CONFIG
# =========================
base_path <- "/scratch/timtd/transcriptomes/MMETSP/fasta/edger_anlysis"

# Input objects (created earlier)
targets <- list(
  diatom = file.path(base_path, "edger_dge_diatom.rds")
)


# Output directory
outdir <- file.path(base_path, "DE_edgeR_QL")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Desired condition order (must match y$samples$condition values)
condition_order <- c("init", "control5", "virus5", "control9", "virus9")

# Normalization method
# - "TMM" is standard
# - "TMMwsp" can help if there are many zeros (often true for bacteria)
norm_method_by_set <- list(
  bacteria = "TMMwsp",
  diatom   = "TMM"
)

# Manual filtering rule
CPM_MIN <- 1
N_MIN   <- 4   # use what you decided previously

# =========================
# Helper: run edgeR QL pipeline
# =========================
run_edger_ql <- function(y, dataset_name, norm_method) {

  # ---- condition handling + ordering ----
  if (!("condition" %in% colnames(y$samples))) {
    stop(dataset_name, ": y$samples must contain a column named 'condition'")
  }

  y$samples$condition <- as.character(y$samples$condition)

  # keep only the 5 expected conditions (drop anything else)
  y <- y[, y$samples$condition %in% condition_order]
  y <- y[, order(match(y$samples$condition, condition_order))]

  y$samples$condition <- factor(y$samples$condition, levels = condition_order)

  # ---- filtering ----
  keep <- filterByExpr(y, group = y$samples$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # ---- STEP 2: manual CPM filter on top of keep1 ----
   # Compute CPM on y1 (already filtered + norm factors)
   y1 <- calcNormFactors(y, method=norm_method)
cpm_mat <- cpm(y1, normalized.lib.sizes=TRUE)

## keep2 <- rowSums(cpm_mat >= CPM_MIN) >= N_MIN   #dont run manual filtering as QC plots don't improve
#y2 <- y1[keep2,, keep.lib.sizes=FALSE]
#y2 <- calcNormFactors(y2, method=norm_method)  # recompute after subsetting
#logcpm2 <- cpm(y2, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)
  

  # ---- design (no intercept) ----
  design <- model.matrix(~ 0 + condition, data = y1$samples)
  colnames(design) <- gsub("^condition", "", colnames(design))  # nicer names

  # ---- dispersion + QL fit ----
  y3 <- estimateDisp(y1, design)
  fit <- glmQLFit(y3, design, robust = TRUE)

  # ---- contrasts ----
  # Primary comparisons:
  # 1) VIRUS5 vs CONTROL5
  # 2) VIRUS9 vs CONTROL9
  # Optional:
  # 3) CONTROL9 vs CONTROL5 (time in controls)
  # 4) VIRUS9 vs VIRUS5 (progression in virus)
  contrast_matrix <- makeContrasts(
    VIRUS5_vs_CONTROL5 = virus5 - control5,
    VIRUS9_vs_CONTROL9 = virus9 - control9,
    CONTROL9_vs_CONTROL5 = control9 - control5,
    VIRUS9_vs_VIRUS5 = virus9 - virus5,
    levels = design
  )

  # ---- run tests + save outputs ----
  res_list <- list()

  for (cn in colnames(contrast_matrix)) {
    qlf <- glmQLFTest(fit, contrast = contrast_matrix[, cn])

    tt <- topTags(qlf, n = Inf, sort.by = "PValue")$table
    tt$feature_id <- rownames(tt)
    tt <- tt[, c("feature_id", setdiff(colnames(tt), "feature_id"))]

    # Save
    out_csv <- file.path(outdir, paste0(dataset_name, "_", cn, "_edgeR_QL.csv"))
    out_rds <- file.path(outdir, paste0(dataset_name, "_", cn, "_edgeR_QL.rds"))
    write.csv(tt, out_csv, row.names = FALSE)
    saveRDS(list(test = qlf, table = tt), out_rds)

    res_list[[cn]] <- tt
  }

  # ---- save fitted objects + normalized expression for later plotting ----
  # logCPM is convenient for heatmaps/QC later
  logcpm <- cpm(y1, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

  saveRDS(y1,   file.path(outdir, paste0(dataset_name, "_DGE_filtered_norm.rds")))
  saveRDS(fit, file.path(outdir, paste0(dataset_name, "_QLfit.rds")))
  saveRDS(logcpm, file.path(outdir, paste0(dataset_name, "_logCPM.rds")))

  # ---- compact QC PDF (post-filter+norm) ----
  qc_pdf <- file.path(outdir, paste0(dataset_name, "_QC_postFilterNorm.pdf"))
  pdf(qc_pdf, width = 10, height = 8)

  par(mfrow = c(2, 2), mar = c(8, 4, 4, 1))

  # library sizes + effective library sizes
  barplot(y1$samples$lib.size / 1e6, las = 2,
          names.arg = y1$samples$condition,
          main = paste0(dataset_name, ": raw lib.size (M)"),
          ylab = "Millions")

  eff_lib <- (y1$samples$lib.size * y1$samples$norm.factors) / 1e6
  barplot(eff_lib, las = 2,
          names.arg = y1$samples$condition,
          main = paste0(dataset_name, ": effective lib.size (M)"),
          ylab = "Millions")

  plotMDS(y1, labels = y1$samples$condition,
          col = as.numeric(y1$samples$condition),
          main = paste0(dataset_name, ": MDS (filtered + ", norm_method, ")"))

  plotBCV(y3, main = paste0(dataset_name, ": BCV (common=", signif(y3$common.dispersion, 3), ")"))

  dev.off()

  message(dataset_name, ": done. Kept ", nrow(y1), " features after filtering.")
  invisible(res_list)
}

# =========================
# RUN for diatom
# =========================
for (nm in names(targets)) {
  message("\n=== Running edgeR QL DE for: ", nm, " ===")

  y <- readRDS(targets[[nm]])

  norm_method <- norm_method_by_set[[nm]]
  if (is.null(norm_method)) norm_method <- "TMM"

  run_edger_ql(y, dataset_name = nm, norm_method = norm_method)
}

message("\nAll done. Results in: ", outdir)


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

####################################
###### ANNOTATION AND ANALYSIS #####
####################################


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})


de_path   <- file.path(base_path, "DE_edgeR_QL")

# Choose dataset: "diatom"
DATASET <- "diatom"

# Load DE result
res <- read.csv(file.path(de_path,
          paste0(DATASET, "_VIRUS5_vs_CONTROL5_edgeR_QL.csv")))

# Load annotation
annot <- readRDS(file.path(base_path, "eggnog_full_named_matched.rds"))

# =========================
# Thresholds
# =========================
FDR_cutoff  <- 0.05
LFC_cutoff  <- 1

# Add regulation category
res <- res %>%
  mutate(
    regulation = case_when(
      FDR < FDR_cutoff & logFC >= LFC_cutoff  ~ "Up",
      FDR < FDR_cutoff & logFC <= -LFC_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )

# =========================
# Volcano plot
# =========================

volcano <- ggplot(res, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color=regulation), alpha=0.6, size=1) +
  scale_color_manual(values=c(
    "Up"="blue",
    "Down"="red",
    "NS"="grey"
  )) +
  geom_vline(xintercept=c(-LFC_cutoff, LFC_cutoff), linetype="dashed") +
  geom_hline(yintercept=-log10(FDR_cutoff), linetype="dashed") +
  theme_minimal() +
  ggtitle(paste0(DATASET, ": VIRUS9 vs CONTROL9")) +
  xlab("log2 Fold Change") +
  ylab("-log10(FDR)")

ggsave(file.path(de_path,
       paste0(DATASET, "_VIRUS5_vs_CONTROL5_volcano.pdf")),
       volcano, width=8, height=6)
	   
# =========================
# Annotate DE features
# =========================

res_annot <- res %>%
  left_join(annot, by=c("feature_id"="X.query"))

# Save full annotated table
write.csv(res_annot,
          file.path(de_path,
            paste0(DATASET, "_VIRUS5_vs_CONTROL5_annotated.csv")),
          row.names=FALSE)

# Save only significant features
sig_res <- res_annot %>%
  filter(FDR < FDR_cutoff & abs(logFC) >= LFC_cutoff)

write.csv(sig_res,
          file.path(de_path,
            paste0(DATASET, "_VIRUS5_vs_CONTROL5_SIGNIFICANT.csv")),
          row.names=FALSE)

message("Volcano + annotation completed for ", DATASET)
