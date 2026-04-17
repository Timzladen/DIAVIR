# Load required packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(Biostrings); packageVersion("Biostrings")


# Read input files
seqtab <- readRDS("seqtab.nochim.rds")        # ASV table
taxa <- readRDS("taxa.rds")                  # Taxonomic assignment
metadata <- metadata <- read.delim("C:/Users/TimTD/OneDrive - NIB/Projekti/ARRS-Postdoc/Rezultati/Pn208 progression experiment/MayExp/16S-sequencing/metadata.txt")  # Sample metadata

# Create phyloseq components
OTU <- otu_table(seqtab, taxa_are_rows = FALSE)
TAX <- tax_table(taxa)
SAMPLES <- sample_data(metadata)
rownames(SAMPLES)<- SAMPLES$SampleID

ps<-phyloseq(otu_table(OTU), tax_table(TAX), sample_data(SAMPLES))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
# Replace NA Genus labels in the taxonomy table before agglomeration
tax_table(ps)[, "Genus"] <- as.character(tax_table(ps)[, "Genus"])
tax_table(ps)[is.na(tax_table(ps)[, "Genus"]), "Genus"] <- "Unclassified"
tax_table(ps)[, "Class"] <- as.character(tax_table(ps)[, "Class"])
tax_table(ps)[is.na(tax_table(ps)[, "Class"]), "Class"] <- "Unclassified"





# Ensure taxonomy columns are character
tax_table(ps)[, "Phylum"] <- as.character(tax_table(ps)[, "Phylum"])
tax_table(ps)[is.na(tax_table(ps)[, "Phylum"]), "Phylum"] <- "Unclassified"

# 🧹 Step 1: REMOVE Cyanobacteria (chloroplast-assigned 16S)
ps_filtered <- subset_taxa(ps, Phylum != "Cyanobacteria" & Phylum != "Cyanobacteriia")


##############################
### Specific taxon analyses   
##############################

#rhodobacteriacea fasta retrieval
rhodobac<- subset_taxa(ps_filtered, Family == "Rhodobacteraceae")
writeXStringSet(refseq(rhodobac), "rhodobacteriacea.fasta")

rhodobac_tax<- as.data.frame(tax_table(rhodobac))
rhodobac_tax$Species<- make.unique(as.character(rhodobac_tax$Species))
tax_table(rhodobac)<-tax_table(rhodobac_tax)
# Step 2: Agglomerate first
ps_class <- tax_glom(ps_filtered, taxrank = "Class")
ps_genus <- tax_glom(ps_filtered, taxrank = "Genus")

# Step 3: THEN normalize
ps_rel_class <- transform_sample_counts(ps_class, function(x) x / sum(x))
ps_rel_genus <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Step 4: Melt to data frame
rel_df_class <- psmelt(ps_rel_class)
rel_df_genus <- psmelt(ps_rel_genus)

# ---- GENUS: Summarize per sample, then compute means + SD ----
genus_sample_sum <- rel_df_genus %>%
  group_by(Sample, Genus, Treatment, Day) %>%
  summarise(sample_abund = sum(Abundance), .groups = "drop")

summary_day_df_genus <- genus_sample_sum %>%
  group_by(Genus, Day, Treatment) %>%
  summarise(
    mean_rel_abund = mean(sample_abund),
    sd_rel_abund = sd(sample_abund),
    .groups = "drop"
  )

# ---- CLASS: Summarize per sample, then compute means + SD ----
class_sample_sum <- rel_df_class %>%
  group_by(Sample, Class, Treatment, Day) %>%
  summarise(sample_abund = sum(Abundance), .groups = "drop")

summary_day_df_class <- class_sample_sum %>%
  group_by(Class, Day, Treatment) %>%
  summarise(
    mean_rel_abund = mean(sample_abund),
    sd_rel_abund = sd(sample_abund),
    .groups = "drop"
  )
# Separate initial values
# Check what your initial rows look like
initial_df_class <- summary_day_df_class %>%
  filter(Treatment == "Initial") %>%
  mutate(Day = 0)

initial_df_genus <- summary_day_df_genus %>%
  filter(Treatment == "Initial") %>%
  mutate(Day = 0)

summary_augmented <- summary_day_df_class %>%
  filter(Treatment != "Initial") %>%
  bind_rows(
    initial_df_class %>% mutate(Treatment = "Control"),
    initial_df_class %>% mutate(Treatment = "Virus")
  )

summary_augmented_genus <- summary_day_df_genus %>%
  filter(Treatment != "Initial") %>%
  bind_rows(
    initial_df_genus %>% mutate(Treatment = "Control"),
    initial_df_genus %>% mutate(Treatment = "Virus")
  )

# ---- Identify Top Taxa ----
top_classes <- summary_augmented %>%
  group_by(Class) %>%
  summarise(avg = mean(mean_rel_abund)) %>%
  slice_max(order_by = avg, n = 10) %>%
  pull(Class)

top_genus <- summary_augmented_genus %>%
  group_by(Genus) %>%
  summarise(avg = mean(mean_rel_abund)) %>%
  slice_max(order_by = avg, n = 10) %>%
  pull(Genus)

ltys <- c("Control" = "dashed",
          "Virus"   = "solid")
cols <- c("Control" = "#4DBBD5FF",
          "Virus"   = "#E64B35FF")

# ---- Plot CLASS with error bars ----
pdf("Bacterial_class_RA.pdf", width = 7, height = 6)
ggplot(filter(summary_augmented, Class %in% top_classes),
       aes(x = Day, y = mean_rel_abund, color = Treatment, linetype = Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_rel_abund - sd_rel_abund,
                    ymax = mean_rel_abund + sd_rel_abund),
                width = 0.3, alpha = 0.5) +
  facet_wrap(~ Class, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = cols)+
  scale_linetype_manual(values = ltys) +
  scale_x_continuous(breaks = sort(unique(summary_augmented_genus$Day))) +
  labs(title = "Relative abundance of bacterial classes",
       y = "Mean relative abundance")
dev.off()

# ---- Plot GENUS with error bars ----

pdf("Bacterial_genera_RA.pdf", width = 7, height = 6)
ggplot(filter(summary_augmented_genus, Genus %in% top_genus),
       aes(x = Day, y = mean_rel_abund, color = Treatment, linetype= Treatment,group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_rel_abund - sd_rel_abund,
                    ymax = mean_rel_abund + sd_rel_abund),
                width = 0.3, alpha = 0.5) +
  facet_wrap(~ Genus, scales = "free_y") +
  scale_color_manual(values = cols)+
  scale_linetype_manual(values = ltys) +
  scale_x_continuous(breaks = sort(unique(summary_augmented_genus$Day))) +
  theme_bw() +
  labs(title = "Relative abundance of bacterial genera",
       y = "Mean relative abundance")
dev.off()
# Extract abundance and taxonomy
abundance_df <- psmelt(ps_rel_genus)
abundance_df_class <- psmelt(ps_rel_class)
# Replace NAs in Genus with "Unclassified"
abundance_df$Genus <- as.character(abundance_df$Genus)
abundance_df$Genus[is.na(abundance_df$Genus)] <- "Unclassified"

abundance_df_class$Class <- as.character(abundance_df_class$Class)
abundance_df_class$Class[is.na(abundance_df_class$Class)] <- "Unclassified"

top15_genera <-abundance_df %>%
  group_by(Genus) %>%
  summarise(avg = mean(Abundance)) %>%
  top_n(15, avg) %>%
  pull(Genus)

top15_class <-abundance_df_class %>%
  group_by(Class) %>%
  summarise(avg = mean(Abundance)) %>%
  top_n(15, avg) %>%
  pull(Class)
#Label other genera as "Others"
abundance_df$Genus_grouped <- ifelse(abundance_df$Genus %in% top15_genera,
                                     abundance_df$Genus, "Others")
abundance_df$condition <- sample_data(ps_rel_genus)$Treatment[match(abundance_df$Sample, rownames(sample_data(ps_rel_genus)))]


abundance_df_class$Class_grouped <- ifelse(abundance_df_class$Class %in% top15_class,
                                     abundance_df_class$Class, "Others")
abundance_df_class$condition <- sample_data(ps_rel_class)$Treatment[match(abundance_df_class$Sample, rownames(sample_data(ps_rel_class)))]

# Prepare plotting data
plot_df <- abundance_df %>%
  group_by(Sample, Genus_grouped, condition) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

plot_df_class <- abundance_df_class %>%
  group_by(Sample, Class_grouped, condition) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")


# Define friendly color palette (max 16 colors)
friendly_palette <- c(
  RColorBrewer::brewer.pal(12, "Set3"),
  "#000000", "#999999", "#E69F00", "#56B4E9"
)[1:length(unique(plot_df$Genus_grouped))]


# Desired sample order
desired_order <- c("Int", "D5C1", "D5C2", "D5C3", 
                   "D5V1", "D5V2", "D5V3", 
                   "D9C1", "D9C2", "D9C3", 
                   "D9V1", "D9V2", "D9V3")

# Ensure plot_df$Sample is exactly a factor with the desired levels
plot_df$Sample <- factor(plot_df$Sample, levels = desired_order)	
plot_df_class$Sample <- factor(plot_df_class$Sample, levels = desired_order)	

# Plot
taxa_barplot<-ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Genus_grouped)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_d3(palette="category20c") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      panel.grid.minor.y = element_blank(),
	  legend.key.size = unit(0.4, "cm"),       # smaller legend boxes
      legend.text = element_text(size = 8),    # smaller legend text
      legend.title = element_text(size = 9)    # optional: slightly smaller title
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))+
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Top 15 Genera Across Samples") 
pdf("Genus_barplot.pdf", width = 7, height = 8)
taxa_barplot
dev.off()
taxa_barplot_class<-ggplot(plot_df_class, aes(x = Sample, y = Abundance, fill = Class_grouped)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_d3(palette="category20c") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.minor.y = element_blank(),
        legend.key.size = unit(0.4, "cm"),       # smaller legend boxes
        legend.text = element_text(size = 8),    # smaller legend text
        legend.title = element_text(size = 9)    # optional: slightly smaller title
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))+
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Top 7 Classes Across Samples") 
pdf("Class_barplot.pdf", width = 7, height = 8)
taxa_barplot_class
dev.off()
####------------------ FAMILY LEVELS repeat ----------------
# 2. Agglomerate to Family level
ps_family <- tax_glom(ps_filtered, taxrank = "Family")

# 3. Normalize AFTER agglomeration
ps_rel_family <- transform_sample_counts(ps_family, function(x) x / sum(x))

# 4. Melt into long format
rel_df_family <- psmelt(ps_rel_family)

# ------------------ SUMMARY FOR LINEPLOTS ------------------

# Summarize per-sample
family_sample_sum <- rel_df_family %>%
  group_by(Sample, Family, Treatment, Day) %>%
  summarise(sample_abund = sum(Abundance), .groups = "drop")

# Compute mean & SD per treatment per day
summary_day_df_family <- family_sample_sum %>%
  group_by(Family, Day, Treatment) %>%
  summarise(
    mean_rel_abund = mean(sample_abund),
    sd_rel_abund = sd(sample_abund),
    .groups = "drop"
  )

# ------------------ ADD INITIAL STATE TO BOTH TREATMENTS ------------------

initial_df_family <- summary_day_df_family %>%
  filter(Treatment == "Initial") %>%
  mutate(Day = 0)

summary_augmented_family <- summary_day_df_family %>%
  filter(Treatment != "Initial") %>%
  bind_rows(
    initial_df_family %>% mutate(Treatment = "Control"),
    initial_df_family %>% mutate(Treatment = "Virus")
  )

# ------------------ IDENTIFY TOP FAMILIES ------------------

top_families <- summary_augmented_family %>%
  group_by(Family) %>%
  summarise(avg = mean(mean_rel_abund)) %>%
  slice_max(order_by = avg, n = 10) %>%
  pull(Family)

# ------------------ PLOT TIME SERIES WITH ERROR BARS ------------------

pdf("Bacterial_family_RA_over_time.pdf", width = 7, height = 6)
ggplot(filter(summary_augmented_family, Family %in% top_families),
       aes(x = Day, y = mean_rel_abund, color = Treatment, linetype= Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_rel_abund - sd_rel_abund,
                    ymax = mean_rel_abund + sd_rel_abund),
                width = 0.3, alpha = 0.5) +
  facet_wrap(~ Family, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = cols)+
  scale_linetype_manual(values = ltys) +
  scale_x_continuous(breaks = sort(unique(summary_augmented_genus$Day))) +
  labs(title = "Relative abundance of bacterial families over time",
       y = "Mean relative abundance", x = "Day")
dev.off()

# ------------------ BARPLOT ACROSS SAMPLES ------------------

# Replace NA Family with "Unclassified"
rel_df_family$Family <- as.character(rel_df_family$Family)
rel_df_family$Family[is.na(rel_df_family$Family)] <- "Unclassified"

# Determine top 15 families for barplot
top15_families <- rel_df_family %>%
  group_by(Family) %>%
  summarise(avg = mean(Abundance)) %>%
  slice_max(order_by = avg, n = 15) %>%
  pull(Family)

# Group others into "Others"
rel_df_family$Family_grouped <- ifelse(rel_df_family$Family %in% top15_families,
                                       rel_df_family$Family, "Others")

# Assign condition label
rel_df_family$condition <- sample_data(ps_rel_family)$Treatment[match(rel_df_family$Sample, rownames(sample_data(ps_rel_family)))]

# Summarize for stacked barplot
plot_df_family <- rel_df_family %>%
  group_by(Sample, Family_grouped, condition) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Reorder samples (if needed)
desired_order <- c("Int", "D5C1", "D5C2", "D5C3", 
                   "D5V1", "D5V2", "D5V3", 
                   "D9C1", "D9C2", "D9C3", 
                   "D9V1", "D9V2", "D9V3")

plot_df_family$Sample <- factor(plot_df_family$Sample, levels = desired_order)

# Plot barplot

cols <- c("Control" = "#4DBBD5FF",
          "Virus"   = "#E64B35FF")
pdf("Bacterial_family_barplot.pdf", width = 7, height = 8)
ggplot(plot_df_family, aes(x = Sample, y = Abundance, fill = Family_grouped)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_d3(palette = "category20c") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.background = element_rect(fill = "white", colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  ylab("Relative Abundance") +
  xlab("Sample") +
  ggtitle("Top 15 Bacterial Families Across Samples")
dev.off()



###################################
###    Diversity analysis 
###################################


alpha_df <- estimate_richness(ps, measures = c("Shannon", "Observed"))
alpha_df$condition <- sample_data(ps)$condition  

# Plot alpha diversity
ggplot(alpha_df, aes(x = condition, y = Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_minimal() +
  ggtitle("Shannon Diversity by Group")
  
###t-test shannon
alpha_df$type<- c(rep("control", 3), rep("virus", 3), rep("control", 3), rep("virus", 3), "init")
alpha_df_noinit<- alpha_df[-13,]

# Kruskal-Wallis test
kruskal_shannon <- kruskal.test(Shannon ~ condition, data = alpha_df)
print(kruskal_shannon)



t_shannon <- t.test(Shannon ~ type, data = alpha_df_noinit)
t_observed <- t.test(Observed ~ type, data = alpha_df_noinit)
print(t_shannon)
print(t_observed)



#        Welch Two Sample t-test
#
#data:  Shannon by type
#t = 2.733, df = 8.9955, p-value = 0.02312
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# 0.05646811 0.59928346
#sample estimates:
#mean in group control   mean in group virus
#             1.564858              1.236982
#
#
#        Welch Two Sample t-test
#
#data:  Observed by type
#t = 1.57, df = 8.7603, p-value = 0.1518
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -2.45791 13.45791
#sample estimates:
#mean in group control   mean in group virus
#             59.16667              53.66667

# Beta diversity (Bray-Curtis)
ord <- ordinate(ps, method = "PCoA", distance = "bray")
p1 <- plot_ordination(ps, ord, color = "Group") +
  geom_point(size = 3) +
  ggtitle("PCoA - Bray-Curtis") +
  theme_minimal()
print(p1)

################################
###   PCoA and ANOSIM 
################################


# ASV matrix from phyloseq object
asv_mat <- as(otu_table(ps), "matrix")

# Make sure samples are rows
if (taxa_are_rows(ps)) {
  asv_mat <- t(asv_mat)
}

# Metadata
meta_df <- data.frame(sample_data(ps))
meta_df$SampleID <- rownames(meta_df)
meta_df <- meta_df[match(rownames(asv_mat), rownames(meta_df)), , drop = FALSE]

stopifnot(all(rownames(asv_mat) == rownames(meta_df)))

# Combined group variable
meta_df$Group <- with(meta_df, ifelse(Treatment == "Initial",
                                      "Initial",
                                      paste(Treatment, Day, sep = "_")))

meta_df$Group <- factor(meta_df$Group,
                        levels = c("Initial", "Control_5", "Virus_5", "Control_9", "Virus_9"))

# Bray-Curtis on raw ASV counts
bray_dist <- vegdist(asv_mat, method = "bray")

# PCoA
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)
eig_vals <- pcoa_res$eig
var_exp <- eig_vals / sum(eig_vals[eig_vals > 0]) * 100

pcoa_df <- data.frame(
  SampleID = rownames(asv_mat),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2]
) %>%
  left_join(meta_df, by = "SampleID")

p <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = Treatment, shape = as.factor(Day))) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCoA based on ASV abundance",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)"),
    shape = "Day"
  )

print(p)

anosim_res <- anosim(bray_dist, grouping = meta_df$Group, permutations = 9999)
print(anosim_res)

pairwise_anosim <- function(dist_mat, grouping, permutations = 9999, p.adjust.method = "BH") {
  grouping <- as.factor(grouping)
  levs <- levels(grouping)
  combs <- combn(levs, 2, simplify = FALSE)
  
  results <- lapply(combs, function(pair) {
    keep <- grouping %in% pair
    
    sub_dist <- as.dist(as.matrix(dist_mat)[keep, keep])
    sub_group <- droplevels(grouping[keep])
    
    fit <- anosim(sub_dist, grouping = sub_group, permutations = permutations)
    
    data.frame(
      group1 = pair[1],
      group2 = pair[2],
      R = unname(fit$statistic),
      p_value = fit$signif,
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results)
  results_df$p_adj <- p.adjust(results_df$p_value, method = p.adjust.method)
  results_df <- results_df[order(results_df$p_value), ]
  rownames(results_df) <- NULL
  results_df
}

# ---------- Run pairwise ANOSIM ----------
pairwise_res <- pairwise_anosim(
  dist_mat = bray_dist,
  grouping = meta_df$Group,
  permutations = 9999,
  p.adjust.method = "BH"
)

####################
###     DESEQ
####################
library(phyloseq)
library(DESeq2)
ps_sub2 <- subset_samples(ps_genus, condition %in% c("initial", "control5"))

dds2 <- phyloseq_to_deseq2(ps_sub2, ~ condition)
dds2 <- DESeq(dds2)
res2 <- results(dds2)
summary(res)
# 1. Extract significant taxa from DESeq2 results
resLFC2 <- lfcShrink(dds2, coef = "condition_initial_vs_control5", type = "apeglm")
sig_res2 <- resLFC2[which(resLFC2$padj < 0.1 & !is.na(resLFC2$padj)), ]  # or use a different threshold
sig_taxa2 <- rownames(sig_res2)

# 2. Extract taxonomy table from phyloseq object
tax_table_df2 <- as.data.frame(tax_table(ps))

# 3. Subset taxonomy to significant taxa
sig_taxonomy2 <- tax_table_df[sig_taxa2, ]

# 4. Combine results with taxonomy for export or inspection
sig_combined2 <- cbind(as.data.frame(sig_res2), sig_taxonomy2)
sig_combined2_nohead<- sig_combined2
rownames(sig_combined2_nohead)<- NULL
# 5. (Optional) View or save
sig_combined2_nohead
write.csv(sig_combined2_nohead, "DE_taxa_init_cont5.csv")


library(dplyr)

# 1. Filter the taxa of interest (Polaribacter or Flavobacteriaceae)
plot_df <- sig_combined %>%
  filter(Genus == "Polaribacter" | Family == "Rhodobacteraceae") %>%
  mutate(TaxonLabel = ifelse(!is.na(Genus), Genus, Family))

# 2. Plot log2 fold changes
ggplot(plot_df, aes(x = reorder(TaxonLabel, log2FoldChange), 
                    y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("red", "steelblue"),
                    labels = c("Down in virus", "Up in virus"),
                    name = "Direction") +
  theme_minimal() +
  labs(title = "Significantly Altered Taxa (Virus vs Control)",
       x = "Taxon",
       y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
