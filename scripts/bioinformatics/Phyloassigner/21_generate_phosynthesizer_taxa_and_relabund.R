library(dplyr)

# ============================================================
# M197 V1-V2
# CYANOBACTERIA + PLASTIDS
# AND EUKARYOTIC PLASTIDS ONLY
#
# Uses PhyloAssigner BEST placements
#
# Produces 4 files:
#   1. CyanoPlastid_relative_abundance.csv
#   2. CyanoPlastid_taxonomy.csv
#   3. Plastid_relative_abundance.csv
#   4. Plastid_taxonomy.csv
#
# Relative abundance:
#   Rows = sequencing samples
#   Columns = ASV sequences
#   Each sample sums to 1 (when group reads > 0)
# ============================================================


# ============================================================
# 1. LOAD DATA
# ============================================================

seqtab <- readRDS("seqtab.nochim.rds")

asv_map <- read.delim(
  "M197_V1V2_ASV_mapping.tsv",
  stringsAsFactors = FALSE
)

tax <- read.csv(
  "07_PhyloAssigner/M197_V1V2_PhyloAssigner_master_taxonomy.csv",
  stringsAsFactors = FALSE
)

# Rename ASV ID to match mapping table
tax <- tax %>%
  rename(ASV = ASV_id)


# ============================================================
# 2. DEFINE TAXONOMIC GROUPS
# ============================================================

cyanobacteria_groups <- c("cya")

plastid_groups <- c(
  "vir",
  "pry",
  "str",
  "rap",
  "cry",
  "rho",
  "din",
  "pav",
  "gla",
  "new_euk_A",
  "new_euk_B",
  "new_euk_C"
)


# ============================================================
# 3. PREPARE TAXONOMY
# ============================================================

tax <- tax %>%
  mutate(
    Broad_group = sub(
      ";.*", "", cyanoplastid_best
    ),

    Group = case_when(
      Broad_group %in% cyanobacteria_groups ~
        "Cyanobacteria",

      Broad_group %in% plastid_groups ~
        "Plastid",

      TRUE ~ NA_character_
    )
  ) %>%
  left_join(
    asv_map %>% select(ASV, Sequence),
    by = "ASV"
  )

stopifnot(!anyDuplicated(tax$ASV))
stopifnot(!anyDuplicated(asv_map$Sequence))
stopifnot(!anyNA(asv_map$ASV))


# ============================================================
# 4. FUNCTION TO GENERATE MATCHED FILES
# ============================================================

make_tables <- function(tax_subset, prefix) {

  # ----------------------------------------
  # Taxonomy table
  # ----------------------------------------

  taxonomy <- tax_subset %>%
    select(
      Sequence,
      Group,
      Broad_group,

      # Primary best placements
      Cyanoplastid_Best = cyanoplastid_best,
      Cyanobacteria_Best = cyanobacteria_best,
      Stramenopiles_Best = stramenopiles_best,
      Pelagophyceae_Best = pelagophyceae_best,
      Dictyochophyceae_Best = dictyochophyceae_best,
      Viridiplantae_Best = viridiplantae_best,

      # Most specialized available best placement
      Final_Best = final_best_taxonomy,
      Final_Database = final_best_database,

      # LCA assignments retained for reference
      Cyanoplastid_LCA = cyanoplastid_LCA,
      Cyanobacteria_LCA = cyanobacteria_LCA,
      Stramenopiles_LCA = stramenopiles_LCA,
      Pelagophyceae_LCA = pelagophyceae_LCA,
      Dictyochophyceae_LCA = dictyochophyceae_LCA,
      Viridiplantae_LCA = viridiplantae_LCA,
      Final_LCA = final_LCA_taxonomy
    )

  if (anyNA(taxonomy$Sequence)) {
    stop("Some classified ASVs are missing sequences.")
  }

  if (anyDuplicated(taxonomy$Sequence)) {
    stop("Duplicate sequences in taxonomy table.")
  }

  # ----------------------------------------
  # Relative abundance table
  # ----------------------------------------

  # Select sequences belonging to this group
  selected_sequences <- taxonomy$Sequence

  # Confirm all sequences occur in DADA2 table
  missing <- setdiff(
    selected_sequences,
    colnames(seqtab)
  )

  if (length(missing) > 0) {
