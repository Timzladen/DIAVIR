#!/usr/bin/env Rscript

# =========================
# 02_edgeR_QC.R
# edgeR QC and normalization diagnostics
# =========================

library(edgeR)

# ---- choose object ----
base_path <- "/scratch/timtd/transcriptomes/MMETSP/fasta/edger_anlysis"

# Options:
# "edger_dge_all.rds"


OBJECT_TO_LOAD <- "edger_dge_all.rds"

y <- readRDS(file.path(base_path, OBJECT_TO_LOAD))
condition_order <- c("init", "control5", "virus5", "control9", "virus9")
# Reorder samples by condition order
ord <- order(match(y$samples$condition, condition_order))
y <- y[, ord]
# Set condition as ordered factor
y$samples$condition <- factor(y$samples$condition, levels = condition_order)

# Output file
pdf(file.path(base_path, paste0("QC_", gsub(".rds","",OBJECT_TO_LOAD), ".pdf")),
    width = 10, height = 16)

par(mfrow=c(4,2))

# =========================
# 1) Raw library sizes
# =========================
barplot(y$samples$lib.size / 1e6,
        las=2,
		names.arg = y$samples$condition,
        ylab="Library size (millions)",
        main="Raw library sizes")

# =========================
# 2) LogCPM BEFORE normalization
# =========================
logcpm_raw <- cpm(y, log=TRUE, prior.count=1)

boxplot(logcpm_raw,
        las=2,
        outline=FALSE,
		names.arg = y$samples$condition,
        main="LogCPM (raw library sizes)",
        ylab="log2 CPM")

# =========================
# 3) TMM normalization
# =========================
y <- calcNormFactors(y, method="TMM")

barplot(y$samples$norm.factors,
        las=2,
		names.arg = y$samples$condition,
        main="TMM normalization factors",
        ylab="Norm factor")

# =========================
# 4) LogCPM AFTER normalization
# =========================
logcpm_norm <- cpm(y, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)

boxplot(logcpm_norm,
        las=2,
		names.arg = y$samples$condition,
        outline=FALSE,
        main="LogCPM (after TMM normalization)",
        ylab="log2 CPM")
		
# =========================
# 5) Boxplots after filtration
# =========================
keep <- filterByExpr(y, group=y$samples$condition)
y_filt <- y[keep,, keep.lib.sizes=FALSE]
# logCPM BEFORE normalization (filtered)
logcpm_before <- cpm(y_filt, log=TRUE, prior.count=1, normalized.lib.sizes=FALSE)

# apply TMM
y_f <- calcNormFactors(y_filt, method="TMM")

# logCPM AFTER normalization (filtered)
logcpm_after <- cpm(y_filt, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)


boxplot(logcpm_before, las=2, outline=FALSE, names.arg = y$samples$condition,
        main="Diatoms (filtered): logCPM pre-TMM", ylab="log2 CPM")

boxplot(logcpm_after, las=2, outline=FALSE,names.arg = y$samples$condition,
        main="Diatoms (filtered): logCPM post-TMM", ylab="log2 CPM")

barplot(y_filt$samples$norm.factors, las=2,names.arg = y$samples$condition,
        main="TMM norm factors (filtered)", ylab="Norm factor")




# =========================
# 6) MDS plot (separate page)
# =========================

plotMDS(y_f,
        labels=colnames(y_f),
        col=as.numeric(y_f$samples$condition),
		main="MDS (filtered + TMM)")

legend("topright",
       legend=levels(y_f$samples$condition),
       col=1:length(levels(y_f$samples$condition)),
       pch=16)

dev.off()		




# =========================
# 5) MDS plot (separate page)
# =========================
pdf(file.path(base_path, paste0("QC_MDS_", gsub(".rds","",OBJECT_TO_LOAD), ".pdf")),
    width=7, height=7)

plotMDS(y,
        labels=colnames(y),
        col=as.numeric(as.factor(y$samples$condition)))

dev.off()

# =========================
# 6) Mean-variance (BCV) trend
# =========================
# Filtering lowly expressed genes first is recommended for this

y_filt <- calcNormFactors(y_filt)

y_filt <- estimateDisp(y_filt, design=model.matrix(~ condition, data=y_filt$samples))

pdf(file.path(base_path, paste0("QC_BCV_", gsub(".rds","",OBJECT_TO_LOAD), ".pdf")),
    width=7, height=7)

plotBCV(y_filt)

dev.off()

message("QC completed for: ", OBJECT_TO_LOAD)

#### density plots

DATASET <- "all"  # "diatom" also works

rds_path <- file.path(base_path, paste0("edger_dge_", DATASET, ".rds"))

condition_order <- c("init", "control5", "virus5", "control9", "virus9")

# Manual filter parameters (start here, adjust if needed)
# Rule: keep genes with CPM >= CPM_MIN in at least N samples
# For bacteria, try CPM_MIN=1 and N=3 or N=4; for diatoms maybe N=3.
CPM_MIN <- 1
N_MIN   <- 3

# Normalization method
norm_method <- if (DATASET == "all") "TMM" else "TMMwsp"

# ---- load and order ----
y <- readRDS(rds_path)
y$samples$condition <- factor(as.character(y$samples$condition), levels=condition_order)
y <- y[, order(match(y$samples$condition, condition_order))]

# ---- STEP 0: logCPM pre-filter, pre-normalization ----
logcpm0 <- cpm(y, log=TRUE, prior.count=1, normalized.lib.sizes=FALSE)

# ---- STEP 1: filterByExpr (still pre-normalization) ----
keep1 <- filterByExpr(y, group=y$samples$condition)
y1 <- y[keep1,, keep.lib.sizes=FALSE]

# normalize AFTER filtering
y1 <- calcNormFactors(y1, method=norm_method)
logcpm1 <- cpm(y1, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)

# ---- STEP 2: manual CPM filter on top of keep1 ----
# Compute CPM on y1 (already filtered + norm factors)
cpm_mat <- cpm(y1, normalized.lib.sizes=TRUE)

keep2 <- rowSums(cpm_mat >= CPM_MIN) >= N_MIN
y2 <- y1[keep2,, keep.lib.sizes=FALSE]
y2 <- calcNormFactors(y2, method=norm_method)  # recompute after subsetting
logcpm2 <- cpm(y2, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)

message(DATASET, " filtering summary:")
message("  features pre-filter:           ", nrow(y))
message("  after filterByExpr:            ", nrow(y1))
message("  after filterByExpr + manual:   ", nrow(y2))
message("  manual rule: CPM >= ", CPM_MIN, " in >= ", N_MIN, " samples")

# ---- density plotting ----
cond_cols <- setNames(seq_along(condition_order), condition_order)
sample_cols0 <- cond_cols[as.character(y$samples$condition)]
sample_cols1 <- cond_cols[as.character(y1$samples$condition)]
sample_cols2 <- cond_cols[as.character(y2$samples$condition)]

dens_panel <- function(logcpm, cols, title) {
  dl <- apply(logcpm, 2, density, na.rm=TRUE)
  xlim <- range(sapply(dl, `[[`, "x"))
  ylim <- range(sapply(dl, `[[`, "y"))
  plot(dl[[1]], col=cols[1], lwd=1.5, xlim=xlim, ylim=ylim,
       main=title, xlab="log2 CPM", ylab="Density")
  if (length(dl) > 1) for (i in 2:length(dl)) lines(dl[[i]], col=cols[i], lwd=1.5)
  legend("topright", legend=levels(y$samples$condition),
         col=cond_cols[levels(y$samples$condition)], lwd=2, bty="n")
}

out_pdf <- file.path(base_path, paste0("QC_density_filter_tuning_", DATASET, ".pdf"))
pdf(out_pdf, width=14, height=5)
par(mfrow=c(1,3), mar=c(4,4,3,1))

dens_panel(logcpm0, sample_cols0, "Pre-filter, no normalization")
dens_panel(logcpm1, sample_cols1, paste0("After filterByExpr + ", norm_method))
dens_panel(logcpm2, sample_cols2, paste0("After + manual CPM filter (", CPM_MIN, " in ", N_MIN, ")"))

dev.off()

message("Wrote: ", out_pdf)



# ---- Plot ----
pdf(file.path(base_path, paste0("QC_boxplots_manual_", DATASET, ".pdf")),
    width=12, height=6)

par(mfrow=c(1,2), mar=c(10,4,4,1))

boxplot(logcpm0,
        names=y$samples$condition,
        las=2,
        outline=FALSE,
        main="Manual filter: pre-normalization",
        ylab="log2 CPM")

boxplot(logcpm2,
        names=y2$samples$condition,
        las=2,
        outline=FALSE,
        main=paste0("Manual filter + ", norm_method),
        ylab="log2 CPM")

dev.off()
