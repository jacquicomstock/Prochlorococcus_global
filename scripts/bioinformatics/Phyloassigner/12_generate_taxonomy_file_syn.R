library(dplyr)

# ASV ID <-> sequence mapping
asv_map <- read.delim(
  "M197_V1V2_ASV_mapping.tsv",
  stringsAsFactors = FALSE
)

# Cyanobacteria LCA assignments
cyano_lca <- read.delim(
  "07_PhyloAssigner/03_cyanobacteria/LCA_placements.tsv",
  stringsAsFactors = FALSE
)

# Cyanobacteria best placements
cyano_best <- read.delim(
  "07_PhyloAssigner/03_cyanobacteria/best_placements.tsv",
  stringsAsFactors = FALSE
)

# Select Synechococcus ASVs based on conservative LCA placement
syn_lca <- cyano_lca %>%
  filter(grepl("MarSyn|Syn5_3", taxon)) %>%
  select(
    ASV = ASV_id,
    Syn_LCA = taxon
  )

# Get corresponding best placements
syn_best <- cyano_best %>%
  filter(ASV_id %in% syn_lca$ASV) %>%
  select(
    ASV = ASV_id,
    Syn_best = taxon,
    likelihood_weight = like_weight_ratio
  )

# Combine with sequences
syn_assignments <- syn_lca %>%
  left_join(syn_best, by = "ASV") %>%
  left_join(asv_map, by = "ASV") %>%
  select(
    ASV,
    Sequence,
    Syn_LCA,
    Syn_best,
    likelihood_weight
  )

# QC
nrow(syn_assignments)
table(syn_assignments$Syn_LCA)

# Save
write.csv(
  syn_assignments,
  "M197_V1V2_Synechococcus_PhyloAssigner.csv",
  row.names = FALSE
)
