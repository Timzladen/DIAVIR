### edgeR diatom plotting ####
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(forcats)

df <- read.csv("diatom_VIRUS5_vs_CONTROL5_SIGNIFICANT.csv", stringsAsFactors = FALSE)

# keep a clean annotation field
df <- df %>%
  mutate(
    Description = ifelse(is.na(Description) | Description == "", "-", Description),
    desc_lc = str_to_lower(Description)
  )

# curate functional bins from description
df_fun <- df %>%
  mutate(
    functional_group = case_when(
      str_detect(desc_lc, "chlorophyll a-b binding|fucoxanthin chlorophyll|chlorophyll a c") ~
        "Light-harvesting antenna",

      str_detect(desc_lc, "photosystem|photosynthetic electron transport|cytochrome b6-f|ferredoxin|psbu|manganese-stabilising") ~
        "Photosynthetic core / electron transport",

      str_detect(desc_lc, "chlorophyll biosynthesis|magnesium chelatase|porphobilinogen|uroporphyrinogen|delta-aminolevulinic|geranylgeranyl diphosphate reductase|zeta-carotene isomerase") ~
        "Pigment / tetrapyrrole biosynthesis",

      str_detect(desc_lc, "silicon transporter") ~
        "Silicon transport / frustule",

      str_detect(desc_lc, "na\\+/pi-cotransporter|phosphate transporter|phosphate permease") ~
        "Phosphate transport / P-stress response",

      str_detect(desc_lc, "hsp|heat shock|chaperone|peptidyl-prolyl cis-trans isomerase") ~
        "Proteostasis / heat shock",

      str_detect(desc_lc, "ribosomal|translation|fibrillarin|snrnp|rna recognition motif|elongation|eukaryotic translation initiation factor|histone h3") ~
        "Translation / RNA processing",

      str_detect(desc_lc, "tubulin|microtubule") ~
        "Cytoskeleton",

      str_detect(desc_lc, "isocitrate lyase|malate synthase|gluconeogenesis|glyceraldehyde 3-phosphate dehydrogenase|transketolase|phosphoglycerate kinase|aldolase|acetyl-coa carboxylase|hydroxymethylglutaryl") ~
        "Central carbon / lipid metabolism",

      str_detect(desc_lc, "oxidoreductase|nad binding|ferredoxin") ~
        "Redox metabolism",

      str_detect(desc_lc, "vacuolar|v-?atpase|sec13|sec63|ran-binding|nucleocytoplasmic|tat complex") ~
        "Trafficking / membrane transport",

      Description == "-" ~ "Unknown",

      TRUE ~ "Other annotated"
    )
  )

# summary by group and direction
group_sum <- df_fun %>%
  group_by(functional_group, regulation) %>%
  summarise(
    n_genes = n(),
    median_logFC = median(logFC, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    signed_n = ifelse(regulation == "Down", -n_genes, n_genes)
  )

# order groups by total represented genes
group_order <- df_fun %>%
  count(functional_group, sort = TRUE) %>%
  pull(functional_group)

group_sum <- group_sum %>%
  mutate(functional_group = factor(functional_group, levels = rev(group_order)))

# main plot
p1 <- ggplot(group_sum, aes(x = signed_n, y = functional_group, fill = median_logFC)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey35") +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  scale_x_continuous(labels = abs) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "Median logFC"
  ) +
  labs(
    x = "Number of differentially expressed genes",
    y = NULL,
    title = "Collapsed functional summary of diatom DE genes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p1

ggsave("Diatom_DE_genes_perannotation.pdf", p1, width=8, height=6)


top_labels <- df_fun %>%
  group_by(functional_group, regulation, Description) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(functional_group, regulation) %>%
  slice_max(order_by = n, n = 3, with_ties = FALSE) %>%
  summarise(
    top_annotations = paste0(Description, " (", n, ")", collapse = "; "),
    .groups = "drop"
  )

support_tab <- group_sum %>%
  left_join(top_labels, by = c("functional_group", "regulation"))

support_tab

write.csv(support_tab, "functional_group_summary_diatom_DE.csv", row.names = FALSE)
