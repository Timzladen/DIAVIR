####KEGG/PFAM enrichment
#Rank-based functional enrichment analysis was performed using the fgsea package in R. Genes were ranked according to the signed edgeR quasi-likelihood F statistic (sign(logFC) × F) derived from the differential expression model, ensuring that both effect direction and statistical support were incorporated into the ranking. PFAM and KEGG gene sets were constructed from transcript annotations, and enrichment significance was assessed using the fast preranked GSEA algorithm implemented in fgsea.

#activate conda deseq2_env
library(tidyverse)
library(fgsea)
library(data.table)

df <- read.csv("bacteria_VIRUS5_vs_CONTROL5_annotated.csv")

head(df)


#build ranks using this metric rank=sign(logFC) X -log(Pvalue) - NOT DE MODEL BASED!
df$rank_stat <- sign(df$logFC) * -log10(df$PValue)

ranks <- df$rank_stat
names(ranks) <- df$feature_id
ranks <- sort(ranks, decreasing = TRUE)

#build ranks based on the DE model
ranks <- sign(df$logFC) * df$F
names(ranks) <- df$feature_id
ranks <- sort(ranks, decreasing = TRUE)

df2 <- df %>%
  mutate(
    PFAMs = ifelse(is.na(PFAMs), "", PFAMs),
    has_GGDEF = str_detect(PFAMs, "(^|,)GGDEF(,|$)"),
    has_EAL   = str_detect(PFAMs, "(^|,)EAL(,|$)"),
    has_PilZ  = str_detect(PFAMs, "(^|,)PilZ(,|$)"),
    has_HDGYP = str_detect(PFAMs, "(^|,)HD_GYP(,|$)|(^|,)HD-GYP(,|$)"),
    arch = case_when(
      has_GGDEF & has_EAL ~ "GGDEF_EAL_hybrid",
      has_GGDEF ~ "GGDEF_only",
      has_EAL ~ "EAL_only",
      has_PilZ ~ "PilZ_only",
      has_HDGYP ~ "HDGYP_only",
      TRUE ~ "other"
    )
  )

pathways_arch <- df2 %>%
  filter(arch != "other") %>%
  group_by(arch) %>%
  summarise(genes = list(unique(feature_id)), .groups = "drop") %>%
  { setNames(.$genes, .$arch) }
  
fgsea_arch <- fgsea(
  pathways = pathways_arch,
  stats = ranks,
  minSize = 5,
  maxSize = 500
) %>%
  arrange(padj)
  
head(fgsea_arch)
  
pfam_map <- df %>%
  select(feature_id, PFAMs) %>%
  filter(PFAMs != "-", PFAMs != "") %>%
  separate_rows(PFAMs, sep=",") %>%
  distinct()

pfam_sets <- split(pfam_map$feature_id, pfam_map$PFAMs)
pfam_fgsea <- fgsea(
  pathways = pfam_sets,
  stats = ranks,
  minSize = 10,
  maxSize = 500
)

pfam_fgsea <- pfam_fgsea %>%
  arrange(padj)

head(pfam_fgsea)


kegg_map <- df %>%
  select(feature_id, KEGG_ko) %>%
  filter(KEGG_ko != "-", KEGG_ko != "") %>%
  separate_rows(KEGG_ko, sep=",") %>%
  distinct()
kegg_sets <- split(kegg_map$feature_id, kegg_map$KEGG_ko)
kegg_fgsea <- fgsea(
  pathways = kegg_sets,
  stats = ranks,
  minSize = 3,
  maxSize = 500
)

kegg_fgsea <- kegg_fgsea %>%
  arrange(padj)

head(kegg_fgsea)


# build GO map
go_map <- df %>%
  select(feature_id, GOs) %>%   # change GOs to your real GO column name if needed
  filter(!is.na(GOs), GOs != "-", GOs != "") %>%
  mutate(GOs = str_replace_all(GOs, "\\s+", "")) %>%   # remove spaces
  separate_rows(GOs, sep = ",") %>%
  filter(GOs != "") %>%
  distinct()

# GO gene sets
go_sets <- split(go_map$feature_id, go_map$GOs)

# enrichment
go_fgsea <- fgsea(
  pathways = go_sets,
  stats = ranks,
  minSize = 10,
  maxSize = 500
) %>%
  arrange(padj)

head(go_fgsea)
