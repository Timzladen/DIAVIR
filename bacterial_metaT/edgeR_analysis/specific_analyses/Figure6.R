###################################################################################
### Stacked bar plots of genus composition within DE-relevant COG categories
### using TMM-normalized CPM
###################################################################################

library(edgeR)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readr)

# -----------------------------
# 0. INPUTS
# -----------------------------
base_path <- "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects"

# edgeR object with TMM/TMMwsp already applied
dge <- readRDS(file.path(base_path, "DE_edgeR_QL", "bacteria_DGE_filtered_norm.rds"))

# annotation
annot <- readRDS(file.path(base_path, "eggnog_full_named_matched.rds"))

# sample metadata if needed
samples <- read.csv("/scratch/timtd/transcriptomes/trimmed/salmon_quant/sample_data.csv",
                    stringsAsFactors = FALSE)
cog_desc <- c(
  A = "RNA processing & modification",
  B = "Chromatin structure & dynamics",
  C = "Energy production & conversion",
  D = "Cell cycle control & mitosis",
  E = "Amino acid metabolism & transport",
  F = "Nucleotide metabolism & transport",
  G = "Carbohydrate metabolism & transport",
  H = "Coenzyme metabolism",
  I = "Lipid metabolism",
  J = "Translation & ribosome structure",
  K = "Transcription",
  L = "DNA replication & repair",
  M = "Cell wall/membrane/envelope biogenesis",
  N = "Cell motility",
  O = "Post-translational modification",
  P = "Inorganic ion transport & metabolism",
  Q = "Secondary metabolite biosynthesis",
  R = "General function prediction",
  S = "Function unknown",
  T = "Signal transduction",
  U = "Intracellular trafficking",
  V = "Defense mechanisms",
  W = "Extracellular structures",
  Y = "Nuclear structure",
  Z = "Cytoskeleton"
)
# -----------------------------
# 1. USER-DEFINED: DE COG set
# -----------------------------
de_cog_desc <- c(
  C = "Energy production & conversion",
  D = "Cell cycle control & mitosis",
  E = "Amino acid metabolism & transport",
  F = "Nucleotide metabolism & transport",
  G = "Carbohydrate metabolism & transport",
  H = "Coenzyme metabolism",
  I = "Lipid metabolism",
  J = "Translation & ribosome structure",
  K = "Transcription",
  L = "DNA replication & repair",
  M = "Cell wall/membrane/envelope biogenesis",
  N = "Cell motility",
  O = "Post-translational modification",
  P = "Inorganic ion transport & metabolism",
  Q = "Secondary metabolite biosynthesis",
  T = "Signal transduction",
  U = "Intracellular trafficking",
  V = "Defense mechanisms",
  W = "Extracellular structures",
  Y = "Nuclear structure",
  Z = "Cytoskeleton"
)

# letter -> description mapping must already exist
# e.g. cog_desc <- c(A="RNA processing...", ...)
de_cog_letters <- names(cog_desc)[cog_desc %in% de_cog_desc]

# -----------------------------
# 2. CALCULATE TMM-NORMALIZED CPM
# -----------------------------
# Non-log CPM, normalized by effective library size
cpm_mat <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)

expr_long <- as.data.frame(cpm_mat) %>%
  tibble::rownames_to_column("transcript_id") %>%
  pivot_longer(
    cols = -transcript_id,
    names_to = "Sample",
    values_to = "CPM"
  )

# -----------------------------
# 3. JOIN ANNOTATION + CONDITION
# -----------------------------
plot_df <- expr_long %>%
  left_join(
    annot %>%
      select(transcript_id, genus, COG_category, domain),
    by = "transcript_id"
  ) %>%
  left_join(
    samples %>% select(sample, condition),
    by = c("Sample" = "sample")
  ) %>%
  filter(!is.na(COG_category), !is.na(genus), !is.na(condition)) %>%
  filter(domain == "Bacteria") %>%
  mutate(
    genus = str_trim(genus),
    condition = as.character(condition),
    COG_function = COG_category
  ) %>%
  filter(str_detect(COG_function, "^[A-Z]$")) %>%
  filter(COG_function %in% de_cog_letters)

# -----------------------------
# 4. SUM CPM BY CONDITION × COG × GENUS
# -----------------------------
summary_cpm <- plot_df %>%
  group_by(condition, COG_function, genus) %>%
  summarise(genus_cpm = sum(CPM, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 5. KEEP TOP GENERA, COLLAPSE REST INTO "Other"
# -----------------------------
top15 <- summary_cpm %>%
  group_by(genus) %>%
  summarise(total_cpm = sum(genus_cpm, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_cpm)) %>%
  slice_head(n = 15) %>%
  pull(genus)

summary_cpm_top <- summary_cpm %>%
  mutate(genus = if_else(genus %in% top15, genus, "Other")) %>%
  group_by(condition, COG_function, genus) %>%
  summarise(genus_cpm = sum(genus_cpm, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 6. FACTOR ORDERING
# -----------------------------
condition_order <- c("init", "control5", "control9", "virus5", "virus9")
cog_order <- de_cog_letters


summary_cpm_top <- summary_cpm_top %>%
  mutate(
    condition = factor(condition, levels = condition_order),
    COG_function = factor(COG_function, levels = cog_order)
  ) %>%
  filter(!is.na(condition), !is.na(COG_function))

# optional: order genera by total CPM, OThers on bottom
genus_order <- summary_cpm_top %>%
  group_by(genus) %>%
  summarise(total = sum(genus_cpm, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genus)

# Put "Other" last so it appears at the bottom of the stacked bars
genus_order <- c(setdiff(genus_order, "Other"), "Other")

summary_cpm_top <- summary_cpm_top %>%
  mutate(genus = factor(genus, levels = genus_order))



# -----------------------------
# 7. PLOT: STACKED BARPLOTS
# -----------------------------
library(ggsci)
pdf("DE_edgeR_QL/COG_stackedbars_TMM_CPM.pdf", width = 12, height = 9)

ggplot(summary_cpm_top,
       aes(x = condition, y = genus_cpm, fill = genus)) +
  geom_col(width = 0.8) +
  facet_wrap(~ COG_function,
             ncol = 3,
             scales = "free_y",
             labeller = labeller(COG_function = function(x) unname(cog_desc[x]))) +
  scale_fill_d3(palette = "category20") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Genus composition of selected COG categories (TMM-normalized CPM)",
    y = "Summed TMM-normalized CPM",
    fill = "Genus"
  )

dev.off()
