library(dplyr)

# ASV ID <-> sequence mapping
asv_map <- read.delim(
  "M197_V1V2_ASV_mapping.tsv",
  stringsAsFactors = FALSE
)

# Conservative LCA assignments
pro_lca <- read.delim(
  "07_PhyloAssigner/04_Prochlorococcus/LCA_placements.tsv",
  stringsAsFactors = FALSE
) %>%
  select(
    ASV = ASV_id,
    Pro_LCA = taxon
  )

# Maximum-likelihood best placements
pro_best <- read.delim(
  "07_PhyloAssigner/04_Prochlorococcus/best_placements.tsv",
  stringsAsFactors = FALSE
) %>%
  select(
    ASV = ASV_id,
    Pro_best = taxon,
    likelihood_weight = like_weight_ratio
  )

# Combine
pro_assignments <- asv_map %>%
  inner_join(pro_best, by = "ASV") %>%
  left_join(pro_lca, by = "ASV") %>%
  select(
    ASV,
    Sequence,
    Pro_LCA,
    Pro_best,
    likelihood_weight
  )

# QC
nrow(pro_assignments)
table(pro_assignments$Pro_LCA, useNA = "ifany")

# Save
write.csv(
  pro_assignments,
  "M197_V1V2_Prochlorococcus_PhyloAssigner.csv",
  row.names = FALSE
)
