library(dplyr)
library(tidyr)

# ============================================================
# M197 V1-V2
# Prochlorococcus + Synechococcus clade abundances
# ============================================================

# -----------------------------
# 1. Load data
# -----------------------------

seqtab <- readRDS("seqtab.nochim.rds")

asv_map <- read.delim(
  "M197_V1V2_ASV_mapping.tsv",
  stringsAsFactors = FALSE
)

pro_tax <- read.csv(
  "M197_V1V2_Prochlorococcus_PhyloAssigner.csv",
  stringsAsFactors = FALSE
)

syn_tax <- read.csv(
  "M197_V1V2_Synechococcus_PhyloAssigner.csv",
  stringsAsFactors = FALSE
)

# -----------------------------
# 2. Convert ASV table to long format
# -----------------------------

counts <- as.data.frame(seqtab)

counts$Sample <- rownames(counts)

counts_long <- counts %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Sequence",
    values_to = "Reads"
  ) %>%
  left_join(
    asv_map,
    by = "Sequence"
  )

# Total 16S reads in each sample
sample_totals <- counts_long %>%
  group_by(Sample) %>%
  summarise(
    Total_16S_reads = sum(Reads),
    .groups = "drop"
  )

# ============================================================
# PROCHLOROCOCCUS
# ============================================================

# -----------------------------
# 3. Add Prochlorococcus taxonomy
# -----------------------------

pro_long <- counts_long %>%
  inner_join(
    pro_tax %>% select(ASV, Pro_LCA),
    by = "ASV"
  )

# -----------------------------
# 4. Simplify taxonomy to terminal LCA
# -----------------------------

pro_long <- pro_long %>%
  mutate(
    Pro_clade = sub(
      ".*;",
      "",
      sub(";$", "", Pro_LCA)
    ),
    Pro_clade = ifelse(
      is.na(Pro_LCA),
      "Unresolved",
      Pro_clade
    )
  )

# -----------------------------
# 5. Sum reads by sample and clade
# -----------------------------

pro_clade <- pro_long %>%
  group_by(Sample, Pro_clade) %>%
  summarise(
    Reads = sum(Reads),
    .groups = "drop"
  )

# Total Prochlorococcus reads/sample
pro_totals <- pro_long %>%
  group_by(Sample) %>%
  summarise(
    Total_Pro_reads = sum(Reads),
    .groups = "drop"
  )

# Calculate both relative abundance metrics
pro_clade <- pro_clade %>%
  left_join(pro_totals, by = "Sample") %>%
  left_join(sample_totals, by = "Sample") %>%
  mutate(
    # Fraction of total 16S reads
    RelAbund_16S = Reads / Total_16S_reads,

    # Fraction of all Prochlorococcus reads
    RelAbund_within_Pro = ifelse(
      Total_Pro_reads > 0,
      Reads / Total_Pro_reads,
      NA_real_
    )
  )

# ============================================================
# SYNECHOCOCCUS
# ============================================================

# -----------------------------
# 6. Add Synechococcus taxonomy
# -----------------------------

syn_long <- counts_long %>%
  inner_join(
    syn_tax %>% select(ASV, Syn_LCA),
    by = "ASV"
  )

# -----------------------------
# 7. Simplify taxonomy to terminal LCA
# -----------------------------

syn_long <- syn_long %>%
  mutate(
    Syn_clade = sub(
      ".*;",
      "",
      sub(";$", "", Syn_LCA)
    )
  )

# -----------------------------
# 8. Sum reads by sample and clade
# -----------------------------

syn_clade <- syn_long %>%
  group_by(Sample, Syn_clade) %>%
  summarise(
    Reads = sum(Reads),
    .groups = "drop"
  )

# Total Synechococcus reads/sample
syn_totals <- syn_long %>%
  group_by(Sample) %>%
  summarise(
    Total_Syn_reads = sum(Reads),
    .groups = "drop"
  )

syn_clade <- syn_clade %>%
  left_join(syn_totals, by = "Sample") %>%
  left_join(sample_totals, by = "Sample") %>%
  mutate(
    # Fraction of total 16S reads
    RelAbund_16S = Reads / Total_16S_reads,

    # Fraction of all Synechococcus reads
    RelAbund_within_Syn = ifelse(
      Total_Syn_reads > 0,
      Reads / Total_Syn_reads,
      NA_real_
    )
  )

# ============================================================
# 9. Save results
# ============================================================

write.csv(
  pro_clade,
  "M197_V1V2_Prochlorococcus_clade_abundances.csv",
  row.names = FALSE
)

write.csv(
  syn_clade,
  "M197_V1V2_Synechococcus_clade_abundances.csv",
  row.names = FALSE
)

# Also save ASV-level tables
write.csv(
  pro_long,
  "M197_V1V2_Prochlorococcus_ASV_abundances.csv",
  row.names = FALSE
)

write.csv(
  syn_long,
  "M197_V1V2_Synechococcus_ASV_abundances.csv",
  row.names = FALSE
)

# ============================================================
# 10. QC
# ============================================================

cat("\nProchlorococcus clades:\n")
print(sort(unique(pro_clade$Pro_clade)))

cat("\nSynechococcus clades:\n")
print(sort(unique(syn_clade$Syn_clade)))

cat("\nNumber of Prochlorococcus ASVs:\n")
print(n_distinct(pro_long$ASV))

cat("\nNumber of Synechococcus ASVs:\n")
print(n_distinct(syn_long$ASV))
