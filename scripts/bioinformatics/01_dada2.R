#module load python/miniforge-25.3.0
#conda activate amplicon

library(dada2)

# ============================================================
# 0. PATHS
# ============================================================

# Raw FASTQ files -- READ ONLY / DO NOT MODIFY
raw_path <- "/project/azworden/sequencing_archive/AVITI_AW-9031521_SIDE_A/SO308/00_raw_data"

# All analysis outputs
out_path <- "/project/azworden/users/jcomstock/SO308"

# SILVA 138.2 database directory
silva_dir <- "/project/azworden/databases/SILVA/SILVA138.2"

# Set working directory so intermediate/final output files
# are written to your SO308 analysis directory
setwd(out_path)


# ============================================================
# 1. INVESTIGATE FASTQ QUALITY
# ============================================================

# Identify paired-end FASTQ files
fnFs <- sort(list.files(
  raw_path,
  pattern = "_R1.fastq.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  raw_path,
  pattern = "_R2.fastq.gz$",
  full.names = TRUE
))

# Check how many files were found
length(fnFs)
length(fnRs)

# Extract sample names
sample.names <- sub("_R1.fastq.gz$", "", basename(fnFs))

# Check that forward and reverse files correspond
stopifnot(
  length(fnFs) == length(fnRs),
  length(fnFs) > 0,
  all(
    sub("_R1.fastq.gz$", "", basename(fnFs)) ==
      sub("_R2.fastq.gz$", "", basename(fnRs))
  )
)

# ------------------------------------------------------------
# Aggregate quality profiles
# ------------------------------------------------------------

pdf(
  file.path(out_path, "SO308_forward_quality_profiles.pdf"),
  width = 12,
  height = 8
)

print(plotQualityProfile(fnFs))

dev.off()


pdf(
  file.path(out_path, "SO308_reverse_quality_profiles.pdf"),
  width = 12,
  height = 8
)

print(plotQualityProfile(fnRs))

dev.off()


# ------------------------------------------------------------
# Individual sample quality profiles
# ------------------------------------------------------------

pdf(
  file.path(out_path, "SO308_forward_quality_profiles_by_sample.pdf"),
  width = 12,
  height = 8
)

for (i in seq_along(fnFs)) {
  print(plotQualityProfile(fnFs[i]))
}

dev.off()


pdf(
  file.path(out_path, "SO308_reverse_quality_profiles_by_sample.pdf"),
  width = 12,
  height = 8
)

for (i in seq_along(fnRs)) {
  print(plotQualityProfile(fnRs[i]))
}

dev.off()


# ============================================================
# 2. FILTERED FILE LOCATIONS
# ============================================================

# Filtered FASTQs go in:
# /project/azworden/users/jcomstock/SO308/filtered

filt_path <- file.path(out_path, "filtered")

dir.create(
  filt_path,
  showWarnings = FALSE,
  recursive = TRUE
)

filtFs <- file.path(
  filt_path,
  paste0(sample.names, "_F_filt.fastq.gz")
)

filtRs <- file.path(
  filt_path,
  paste0(sample.names, "_R_filt.fastq.gz")
)


# ============================================================
# 3. TRIM PRIMER REGIONS AND FILTER READS
# ============================================================

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  
  # Remove primer-derived sequence from 5' ends:
  # 27F = 20 nt
  # 338R = 19 nt
  trimLeft = c(20, 19),
  
  # Truncate based on quality profiles
  # Final retained lengths:
  # Forward = 280 - 20 = 260 nt
  # Reverse = 240 - 19 = 221 nt
  truncLen = c(280, 240),
  
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  
  compress = TRUE,
  multithread = TRUE
)

# Examine filtering results
head(out)

# Save filtering statistics
write.table(
  out,
  file.path(out_path, "filtering_stats.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ============================================================
# 4. DEREPLICATE
# ============================================================

derepFs <- derepFastq(
  filtFs,
  verbose = TRUE
)

derepRs <- derepFastq(
  filtRs,
  verbose = TRUE
)

names(derepFs) <- sample.names
names(derepRs) <- sample.names


# ============================================================
# 5. LEARN ERROR RATES
# ============================================================

dadaFs.lrn <- dada(
  derepFs,
  err = NULL,
  selfConsist = TRUE,
  multithread = TRUE
)

errF <- dadaFs.lrn[[1]]$err_out


dadaRs.lrn <- dada(
  derepRs,
  err = NULL,
  selfConsist = TRUE,
  multithread = TRUE
)

errR <- dadaRs.lrn[[1]]$err_out


# Save error models
saveRDS(
  errF,
  file.path(out_path, "errF.rds")
)

saveRDS(
  errR,
  file.path(out_path, "errR.rds")
)


# Save error model plots
pdf(
  file.path(out_path, "SO308_error_rates_forward.pdf"),
  width = 10,
  height = 8
)

print(plotErrors(errF, nominalQ = TRUE))

dev.off()


pdf(
  file.path(out_path, "SO308_error_rates_reverse.pdf"),
  width = 10,
  height = 8
)

print(plotErrors(errR, nominalQ = TRUE))

dev.off()


# ============================================================
# 6. DENOISE
# ============================================================

dadaFs <- dada(
  derepFs,
  err = errF,
  multithread = TRUE
)

dadaRs <- dada(
  derepRs,
  err = errR,
  multithread = TRUE
)


saveRDS(
  dadaFs,
  file.path(out_path, "dadaFs.rds")
)

saveRDS(
  dadaRs,
  file.path(out_path, "dadaRs.rds")
)


# ============================================================
# 7. MERGE PAIRED READS
# ============================================================

mergers <- mergePairs(
  dadaFs,
  derepFs,
  dadaRs,
  derepRs,
  verbose = TRUE
)

saveRDS(
  mergers,
  file.path(out_path, "mergers.rds")
)


# ============================================================
# 8. CREATE SEQUENCE TABLE
# ============================================================

# Remove mock community if there is a sample named exactly "Mock"
mergers.use <- mergers[names(mergers) != "Mock"]

seqtab <- makeSequenceTable(mergers.use)

dim(seqtab)

# Examine merged sequence-length distribution
table(nchar(getSequences(seqtab)))

# Save pre-chimera sequence table
saveRDS(
  seqtab,
  file.path(out_path, "seqtab.rds")
)


# ============================================================
# 9. REMOVE CHIMERAS
# ============================================================

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

saveRDS(
  seqtab.nochim,
  file.path(out_path, "seqtab.nochim.rds")
)


# ============================================================
# 10. TRACK READS THROUGH PIPELINE
# ============================================================

getN <- function(x) sum(getUniques(x))

track <- cbind(
  input = out[, "reads.in"],
  filtered = out[, "reads.out"],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim)
)

rownames(track) <- sample.names

write.table(
  track,
  file.path(out_path, "read_tracking.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

head(track)


# ============================================================
# 11. TAXONOMIC ASSIGNMENT
# ============================================================

taxa <- assignTaxonomy(
  seqtab.nochim,
  silva_train,
  multithread = TRUE
)

saveRDS(
  taxa,
  file.path(out_path, "taxa.rds")
)

# ============================================================
# 12. EXPORT RESULTS
# ============================================================

write.table(
  cbind(t(seqtab.nochim), taxa),
  file.path(out_path, "SO308_seqtab-nochimtaxa.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)

write.table(
  taxa,
  file.path(out_path, "SO308_taxa.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)
