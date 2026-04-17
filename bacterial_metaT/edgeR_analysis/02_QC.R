#!/usr/bin/env Rscript

# =========================
# 02_edgeR_QC.R
# edgeR QC and normalization diagnostics
# =========================

library(edgeR)

# ---- choose object ----
base_path <- "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects"

# Options:
# "edger_dge_all.rds"
# "edger_dge_bacteria.rds"
# "edger_dge_diatom.rds"

OBJECT_TO_LOAD <- "edger_dge_diatom.rds"

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


