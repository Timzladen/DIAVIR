library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidyr)


# -----------------------------
# Input files
# -----------------------------
file_day5 <- "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/bacteria_VIRUS5_vs_CONTROL5_SIGNIFICANT.csv"
file_day9 <- "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/bacteria_VIRUS9_vs_CONTROL9_SIGNIFICANT.csv"
file_control <- "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/bacteria_CONTROL9_vs_CONTROL5_SIGNIFICANT.csv"

# -----------------------------
# Read data
# -----------------------------
day5 <- read_csv(file_day5, show_col_types = FALSE) %>%
  mutate(condition = "5")

day9 <- read_csv(file_day9, show_col_types = FALSE) %>%
  mutate(condition = "9")
  
control<- read_csv(file_control,  show_col_types= FALSE)

# -----------------------------
# Combine and clean
# -----------------------------
df <- bind_rows(day5, day9) %>%
  filter(!is.na(genus), genus != "", !is.na(regulation)) %>%
  mutate(
    regulation = factor(regulation, levels = c("Up", "Down")),
    condition = factor(condition, levels = c("5", "9"))
  )

# -----------------------------
# Keep top 10 genera by total DE genes
# across both conditions combined
# -----------------------------
top_genera <- df %>%
  count(genus, name = "total_genes") %>%
  arrange(desc(total_genes)) %>%
  slice_head(n = 10) %>%
  pull(genus)

plot_df <- df %>%
  filter(genus %in% top_genera) %>%
  count(condition, genus, regulation, name = "n_genes")

# ensure missing Up/Down combinations appear as 0
plot_df <- plot_df %>%
  tidyr::complete(condition, genus, regulation, fill = list(n_genes = 0))

# order genera by total counts within selected top 10
genus_order <- plot_df %>%
  group_by(genus) %>%
  summarise(total = sum(n_genes), .groups = "drop") %>%
  arrange(total) %>%
  pull(genus)

plot_df <- plot_df %>%
  mutate(genus = factor(genus, levels = genus_order))

# -----------------------------
# Plot option 1:
# side-by-side Up/Down bars, facetted by condition
# -----------------------------
p1 <- ggplot(plot_df, aes(x = genus, y = n_genes, fill = regulation)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~ condition, nrow = 1) +
  coord_flip() +
  labs(
    x = "Genus",
    y = "Number of significant genes",
    fill = "Regulation",
    title = "Differentially expressed genes by genus",
    subtitle = "Top 10 genera by total number of DE genes across day 5 and day 9"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    panel.grid.major.y = element_blank()
  )

print(p1)

ggsave(
  filename = "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/DE_genes_by_genus_top10_facet_condition.pdf",
  plot = p1,
  width = 10,
  height = 6
)

# -----------------------------
# Plot option 2:
# stacked bars, facetted by condition
# -----------------------------
p2 <- ggplot(plot_df, aes(x = genus, y = n_genes, fill = regulation)) +
  geom_col(width = 0.75) +
  facet_wrap(~ condition, nrow = 1) +
  coord_flip() +
  labs(
    x = "Genus",
    y = "Number of significant genes",
    fill = "Regulation",
    title = "Differentially expressed genes by genus",
    subtitle = "Stacked view; top 10 genera across day 5 and day 9"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    panel.grid.major.y = element_blank()
  )

print(p2)

ggsave(
  filename = "DE_genes_by_genus_top10_stacked_facet_condition.pdf",
  plot = p2,
  width = 10,
  height = 6
)


###COG plots per Genus
df <- day5 %>%
  filter(!is.na(genus), genus != "", !is.na(regulation)) %>%
  mutate(
    regulation = factor(regulation, levels = c("Up", "Down"))
  )


# -----------------------------
# Genera of interest
# -----------------------------
target_genera <- c(
  "Polaribacter",
  "Fluviicola",
  "Marinobacter",
  "Alteromonas",
  "Muricauda",
  "Crocinitomix"
)

# -----------------------------
# COG category names
# -----------------------------
cog_desc <- c(
  J = "Translation, ribosomal structure and biogenesis",
  A = "RNA processing and modification",
  K = "Transcription",
  L = "Replication, recombination and repair",
  B = "Chromatin structure and dynamics",
  D = "Cell cycle control, cell division, chromosome partitioning",
  Y = "Nuclear structure",
  V = "Defense mechanisms",
  T = "Signal transduction mechanisms",
  M = "Cell wall/membrane/envelope biogenesis",
  N = "Cell motility",
  Z = "Cytoskeleton",
  W = "Extracellular structures",
  U = "Intracellular trafficking, secretion, and vesicular transport",
  O = "Posttranslational modification, protein turnover, chaperones",
  C = "Energy production and conversion",
  G = "Carbohydrate transport and metabolism",
  E = "Amino acid transport and metabolism",
  F = "Nucleotide transport and metabolism",
  H = "Coenzyme transport and metabolism",
  I = "Lipid transport and metabolism",
  P = "Inorganic ion transport and metabolism",
  Q = "Secondary metabolites biosynthesis, transport and catabolism",
  R = "General function prediction only",
  S = "Function unknown",
  X = "Mobilome: prophages, transposons",
  "-" = "Unannotated"
)

# Optional fixed order
cog_order <- c(
  "J","A","K","L","B","D","Y","V","T","M","N","Z","W","U","O",
  "C","G","E","F","H","I","P","Q","R","S","X","-"
)

# -----------------------------
# Clean and summarize
# -----------------------------
plot_df <- df %>%
  filter(
    genus %in% target_genera,
    !is.na(regulation),
    !is.na(COG_category),
    COG_category != ""
  ) %>%
  mutate(
    # Collapse multi-letter COGs to first letter
    COG_category = substr(COG_category, 1, 1),
    COG_category = ifelse(COG_category %in% names(cog_desc), COG_category, "-"),
    regulation = factor(regulation, levels = c("Up", "Down")),
    genus = factor(genus, levels = target_genera)
  ) %>%
  count(genus, COG_category, regulation, name = "n_genes") %>%
  complete(genus, COG_category = cog_order, regulation, fill = list(n_genes = 0)) %>%
  mutate(
    COG_description = paste0(COG_category, ": ", cog_desc[COG_category])
  )

# keep only categories present in at least one direction within a genus
plot_df_nonzero <- plot_df %>%
  group_by(genus, COG_description) %>%
  filter(sum(n_genes) > 0) %>%
  ungroup()

# order COGs globally by total counts
cog_levels <- plot_df_nonzero %>%
  group_by(COG_description) %>%
  summarise(total = sum(n_genes), .groups = "drop") %>%
  arrange(total) %>%
  pull(COG_description)

plot_df_nonzero <- plot_df_nonzero %>%
  mutate(COG_description = factor(COG_description, levels = cog_levels))
# -----------------------------
# Plot 3: faceted by genus and regulation
# -----------------------------
p3 <- ggplot(plot_df, aes(x = COG_description, y = n_genes)) +
  geom_col(width = 0.8) +
  facet_grid(genus ~ regulation, scales = "free_y") +
  coord_flip() +
  labs(
    x = "COG category",
    y = "Number of DE genes",
    title = "Day 5 DE genes by COG category in selected bacterial genera"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8)
  )

print(p3)

ggsave(
  "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/Day5_selected_genera_COG_by_regulation_facetgrid.pdf",
  p1,
  width = 14,
  height = 12
)



# -----------------------------
# Plot 1: stacked Up/Down in same bar
# -----------------------------
# -----------------------------
# Prepare data
# -----------------------------
library(patchwork)

plot_df <- df %>%
  filter(
    genus %in% target_genera,
    !is.na(regulation),
    !is.na(COG_category),
    COG_category != ""
  ) %>%
  mutate(
    COG_category = substr(COG_category, 1, 1),
    COG_category = ifelse(COG_category %in% names(cog_desc), COG_category, "-"),
    regulation = factor(regulation, levels = c("Down", "Up")),  # reversed order
    genus = factor(genus, levels = target_genera)
  ) %>%
  count(genus, COG_category, regulation, name = "n_genes") %>%
  complete(genus, COG_category = cog_order, regulation, fill = list(n_genes = 0))

plot_df_nonzero <- plot_df %>%
  group_by(genus, COG_category) %>%
  filter(sum(n_genes) > 0) %>%
  ungroup()

# global order by total counts, but using letters only
cog_levels <- plot_df_nonzero %>%
  group_by(COG_category) %>%
  summarise(total = sum(n_genes), .groups = "drop") %>%
  arrange(total) %>%
  pull(COG_category)

plot_df_nonzero <- plot_df_nonzero %>%
  mutate(COG_category = factor(COG_category, levels = cog_levels))

p_stacked <- ggplot(plot_df_nonzero, aes(x = COG_category, y = n_genes, fill = regulation)) +
  geom_col(width = 0.8) +
  facet_wrap(~ genus, scales = "free_y", ncol = 2) +
  coord_flip() +
  labs(
    x = "COG category",
    y = "Number of DE genes",
    fill = "Regulation",
    title = "Day 5 DE genes by COG category in selected bacterial genera"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "black"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
	strip.text = element_text(face = "italic"),
    legend.position = "bottom"
  )

# make a text block with the COG descriptions
desc_df <- plot_df_nonzero %>%
  distinct(COG_category) %>%
  mutate(
    COG_category = as.character(COG_category),
    label = paste0(COG_category, " = ", cog_desc[COG_category])
  )
#split into 3 columns  
desc_df <- plot_df_nonzero %>%
  distinct(COG_category) %>%
  mutate(
    COG_category = as.character(COG_category),
    label = paste0(COG_category, " = ", cog_desc[COG_category]),
    idx = row_number(),
    col = (idx - 1) %/% ceiling(n() / 3) + 1,
    row = ave(idx, col, FUN = seq_along)
  ) %>%
  select(-idx)

p_desc <- ggplot(desc_df, aes(x = col, y = -row, label = label)) +
  geom_text(hjust = 0, size = 2.5) +
  xlim(0.7, 3.3) +
  theme_void() +
  theme(
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )
  
  
final_plot <- p_stacked / p_desc + plot_layout(heights = c(4.5, 1.8))

print(final_plot)

ggsave(
  "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/Day5_selected_genera_COG_up_down_d5-d9_stacked.pdf",
  plot = final_plot,
  width = 210,
  height = 297,
  units = "mm"
)
