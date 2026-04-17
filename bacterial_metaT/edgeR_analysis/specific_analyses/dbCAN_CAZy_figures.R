library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# input files. X is the Day (5 or 9)
read_de<- function (x) {
de_file <- file.path("/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL",paste0("bacteria_VIRUS",x,"_vs_CONTROL",x,"_SIGNIFICANT.csv"))
read_csv(de_file, show_col_types = FALSE)
}
dbcan_file <- "/scratch/timtd/transcriptomes/trimmed/assembled/dbcan_output/overview.txt"

de<-read_de(5)
dbcan <- read_tsv(dbcan_file, show_col_types = FALSE)


n0 <- nrow(de)

collapse_empty <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x %in% c("", "-")] <- NA_character_
  x
}

extract_cazy_token <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  
  hits <- stringr::str_extract_all(x, "(AA|CBM|CE|GH|GT|PL)[0-9]+(_[0-9]+)?")
  
  vapply(hits, function(h) {
    h <- unique(h[!is.na(h) & h != ""])
    if (length(h) == 0) NA_character_ else paste(h, collapse = ";")
  }, character(1))
}

dbcan2 <- dbcan %>%
  mutate(
    feature_id = str_remove(`Gene ID`, "\\.p\\d+$"),
    HMMER_raw = if ("HMMER" %in% colnames(.)) collapse_empty(HMMER) else NA_character_,
    DIAMOND_raw = if ("DIAMOND" %in% colnames(.)) collapse_empty(DIAMOND) else NA_character_,
    dbCAN_sub_raw = if ("dbCAN_sub" %in% colnames(.)) collapse_empty(dbCAN_sub) else NA_character_,
    nTools = if ("#ofTools" %in% colnames(.)) suppressWarnings(as.numeric(`#ofTools`)) else NA_real_
  ) %>%
  mutate(
    HMMER_clean = extract_cazy_token(HMMER_raw),
    DIAMOND_clean = extract_cazy_token(DIAMOND_raw),
    dbCAN_sub_clean = extract_cazy_token(dbCAN_sub_raw),
    has_hmmer = !is.na(HMMER_clean),
    has_sub = !is.na(dbCAN_sub_clean),
    has_diamond = !is.na(DIAMOND_clean)
  )

dbcan_best <- dbcan2 %>%
  group_by(feature_id) %>%
  arrange(
    desc(nTools),
    desc(has_hmmer),
    desc(has_sub),
    desc(has_diamond),
    `Gene ID`,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    CAZy_reannot = case_when(
      !is.na(nTools) & nTools >= 2 ~ coalesce(HMMER_clean, dbCAN_sub_clean, DIAMOND_clean),
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    feature_id,
    dbCAN_GeneID = `Gene ID`,
    dbCAN_HMMER = HMMER_raw,
    dbCAN_DIAMOND = DIAMOND_raw,
    dbCAN_sub = dbCAN_sub_raw,
    dbCAN_tools = nTools,
    CAZy_reannot
  )

de_out <- de %>%
  left_join(dbcan_best, by = "feature_id")

stopifnot(nrow(de_out) == n0)

out_file <- "/DATA/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/bacteria_VIRUS9_CONTROL9_SIGNIFICANT_reannotated.csv"
write_csv(de_out, out_file)

cat("Original rows:", n0, "\n")
cat("Output rows:  ", nrow(de_out), "\n")
cat("Wrote:", out_file, "\n")
cat("Rows with new CAZy_reannot:", sum(!is.na(de_out$CAZy_reannot)), "\n")


###PLOTS

library(ggplot2)
library(forcats)

clean_ann <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- str_trim(x)
  x[x %in% c("", "-", "NA", "None")] <- NA_character_
  x
}

extract_valid_cazy <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  
  out <- str_extract_all(x, "(AA|CBM|CE|GH|GT|PL)[0-9]+(_[0-9]+)?")
  out <- lapply(out, unique)
  out
}

#genera with more than 2 DE CAZy genes
target_genera <- c(
  "Polaribacter",
  "Fluviicola",
  "Marinobacter",
  "Alteromonas",
  "Muricauda",
  "Crocinitomix",
  "Arenitalea",
  "Flavobacteriaceae",
  "Roseovarius",
  "Saprospira",
  "Zunongwangia"
)

#----------------------------
# prepare DE direction
#----------------------------
df2 <- de_out %>%
  mutate(
    direction = case_when(
      logFC > 0 ~ "Upregulated",
      logFC < 0 ~ "Downregulated",
      TRUE ~ NA_character_
    ),
    CAZy = clean_ann(CAZy),
    CAZy_reannot = clean_ann(CAZy_reannot)
  ) %>%
  filter(!is.na(direction)) %>%
  #filter(genus %in% target_genera) %>%
  mutate(
    genus = factor(genus),
	#,levels = target_genera),
    has_any_cazy = !is.na(CAZy) | !is.na(CAZy_reannot)
  )

  
#----------------------------
# PLOT 1
# count genes with either CAZy or CAZy_reannot
# if both exist, count gene only once
#----------------------------
plot1_df <- df2 %>%
  filter(has_any_cazy) %>%
  distinct(feature_id, genus, direction) %>%
  count(genus, direction, name = "n_genes")

p1 <- ggplot(plot1_df, aes(x = direction, y = n_genes, fill = direction)) +
  geom_col(width = 0.7) +
  facet_wrap(~ genus, scales = "free_y") +
  labs(
    x = NULL,
    y = "Number of genes",
    title = "Differentially expressed CAZy-annotated genes by genus"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    strip.text = element_text(face = "italic")
  )

print(p1)

# save if wanted. Make sure to change D.
#ggsave("CAZy_gene_counts_up_down_D9.png", p1, width = 5, height = 4, dpi = 600)
ggsave("CAZy_gene_counts_up_down_genus_D5.png", p1, width = 10, height = 8, dpi = 600)

#----------------------------
# PLOT 2
# CAZy family composition among up/down genes
# prefer CAZy_reannot, otherwise use CAZy
#----------------------------
plot2_long <- df2 %>%
  mutate(
    CAZy_final = coalesce(CAZy_reannot, CAZy)
  ) %>%
  filter(!is.na(CAZy_final)) %>%
  mutate(
    cazy_list = extract_valid_cazy(CAZy_final)
  ) %>%
  select(feature_id, direction, cazy_list) %>%
  unnest(cazy_list) %>%
  rename(CAZy_family = cazy_list) %>%
  filter(!is.na(CAZy_family), CAZy_family != "") %>%
  distinct(feature_id, direction, CAZy_family)

# optional: keep top families overall
top_n <- 20

top_fams <- plot2_long %>%
  count(CAZy_family, sort = TRUE) %>%
  #slice_head(n = top_n) %>%
  pull(CAZy_family)

plot2_df <- plot2_long %>%
  filter(CAZy_family %in% top_fams) %>%
  count(direction, CAZy_family, name = "n_genes") %>%
  group_by(CAZy_family) %>%
  mutate(total = sum(n_genes)) %>%
  ungroup() %>%
  mutate(CAZy_family = fct_reorder(CAZy_family, total))

p2 <- ggplot(plot2_df, aes(x = CAZy_family, y = as.integer(n_genes), fill = direction)) +
  geom_col(position = "dodge") +
  coord_flip() +
   scale_y_continuous(breaks = scales::breaks_width(1)) +
  labs(
    x = "CAZy family",
    y = "Number of genes",
    title = paste0("CAZy families in up- and downregulated genes (top ", top_n, ")")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

print(p2)

ggsave("CAZy_family_up_down_D5.png", p2, width = 8, height = 7, dpi = 600)

################################
###   PLOT 3 D5 and D9 combined
################################




library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
file_day5 <- "/DATA/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/bacteria_VIRUS5_CONTROL5_SIGNIFICANT_reannotated.csv"
file_day9 <- "/DATA/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst/edger_objects/DE_edgeR_manualfilt_QL/bacteria_VIRUS9_CONTROL9_SIGNIFICANT_reannotated.csv"

genus_levels <- c(target_genera, "Other")


clean_ann <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- str_trim(x)
  x[x %in% c("", "-", "NA", "None")] <- NA_character_
  x
}

prep_day <- function(file, day_label) {
  read_csv(file, show_col_types = FALSE) %>%
    mutate(
      Day = day_label,
      direction = case_when(
        logFC > 0 ~ "Upregulated",
        logFC < 0 ~ "Downregulated",
        TRUE ~ NA_character_
      ),
      CAZy = clean_ann(CAZy),
      CAZy_reannot = clean_ann(CAZy_reannot),
      has_any_cazy = !is.na(CAZy) | !is.na(CAZy_reannot),
      genus = if_else(genus %in% target_genera, genus, "Other")
    ) %>%
    filter(!is.na(direction), has_any_cazy) %>%
    distinct(feature_id, genus, Day, direction)
}


#----------------------------
# read and combine
#----------------------------
df <- bind_rows(
  prep_day(file_day5, "Day 5"),
  prep_day(file_day9, "Day 9")
) %>%
  mutate(
    genus = factor(genus, levels = genus_levels),
    direction = factor(direction, levels = c("Upregulated", "Downregulated")),
    Day = factor(Day, levels = c("Day 5", "Day 9"))
  )
#----------------------------
# plot data
#----------------------------

plot_df <- df %>%
  count(genus, direction, Day, name = "n_genes")

plot_df <- expand_grid(
  genus = factor(genus_levels, levels = genus_levels),
  direction = factor(c("Upregulated", "Downregulated"),
                     levels = c("Upregulated", "Downregulated")),
  Day = factor(c("Day 5", "Day 9"), levels = c("Day 5", "Day 9"))
) %>%
  left_join(plot_df, by = c("genus", "direction", "Day")) %>%
  mutate(n_genes = replace_na(n_genes, 0))


#----------------------------
# main plot
#----------------------------
p <- ggplot(plot_df, aes(x = genus, y = n_genes, fill = direction)) +
  geom_col(
    aes(group = interaction(direction, Day)),
    position = position_dodge2(width = 0.8, preserve = "single"),
    width = 0.7
  ) +
  facet_wrap(~ Day, nrow = 1, scales = "free_y") +
  labs(
    x = NULL,
    y = "Number of genes",
    title = "CAZy-annotated DE genes by genus",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#----------------------------
# summary table under plot
# total up/down per day
#----------------------------
tab_df <- df %>%
  count(Day, direction, name = "n_genes") %>%
  tidyr::pivot_wider(
    names_from = direction,
    values_from = n_genes,
    values_fill = 0
  ) %>%
  arrange(Day)

table_grob <- tableGrob(
  tab_df,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 11,
    core = list(fg_params = list(hjust = 0.5, x = 0.5)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

title_grob <- textGrob(
  "Total CAZy-annotated DE genes per day",
  gp = gpar(fontsize = 12, fontface = "bold")
)

table_block <- arrangeGrob(
  #title_grob,
  table_grob,
  ncol = 1,
  heights = c(0.2, 1)
)

#----------------------------
# combine plot + table
#----------------------------
final_plot <- plot_grid(
  p,
  table_block,
  ncol = 1,
  rel_heights = c(3.3, 1)
)

print(final_plot)

ggsave(
  "CAZy_genus_day5_day9_with_table.pdf",
  final_plot,
  width = 11,
  height = 8
)
