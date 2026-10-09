# ============================================================
# Create ASV FASTA for PhyloAssigner
# M197 V1-V2
# ============================================================

seqtab <- readRDS("seqtab.nochim.rds")

# ASV sequences are stored as column names
asv_seqs <- colnames(seqtab)

# Assign simple ASV identifiers
asv_ids <- paste0("ASV_", seq_along(asv_seqs))

# Write FASTA
fasta <- as.vector(
  rbind(
    paste0(">", asv_ids),
    asv_seqs
  )
)

writeLines(
  fasta,
  "M197_V1V2_ASVs.fasta"
)

# Save correspondence between ASV IDs and sequences
asv_map <- data.frame(
  ASV = asv_ids,
  Sequence = asv_seqs
)

write.table(
  asv_map,
  "M197_V1V2_ASV_mapping.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
