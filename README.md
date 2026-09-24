# Prochlorococcus_global

## Overview

This repository contains code and associated metadata for analysis of 16S rRNA
gene amplicon sequencing data and flowcytometry data from oceanographic research cruises 
in the following locations: North Pacific, Indian Ocean, Sargasso Sea, South Atlantic, 
Mediterranean, and the South Pacific.

The analyses resolve the global inventory of prochlorococcus throughout the euphotic water column, 
with 16S phylogeny characterized alongside absolute cell counts.

## Sequencing

Samples were sequenced using paired-end Illumina sequencing targeting the
V1–V2 region of the 16S rRNA gene.

Primers:
- Forward: 27F (5'-AGAGTTTGATCMTGGCTCAG-3')
- Reverse: 338R (5'-TGCTGCCTCCCGTAGGAGT-3')

Primer sequences were removed using Cutadapt prior to sequence processing.

## Sequence processing

Amplicon sequence variants (ASVs) were inferred using the DADA2 pipeline in R.

The general workflow was:

1. Remove amplification primers using Cutadapt
2. Inspect read quality
3. Filter and truncate reads
4. Learn forward and reverse error rates
5. Dereplicate reads
6. Infer ASVs using DADA2
7. Merge paired-end reads
8. Remove chimeric sequences
9. Assign taxonomy using [SILVA version]
10. Remove non-target sequences and perform downstream analyses

Forward and reverse reads were truncated to 280 and 240 bp, respectively,
following primer removal.

Additional filtering parameters and software versions are documented in the
analysis scripts.

## Repository structure

    .
    ├── README.md
    ├── LICENSE
    ├── scripts/
    │   ├── 01_cutadapt.sh
    │   ├── 02_dada2.R
    │   └── 03_downstream_analysis.R
    ├── metadata/
    │   └── sample_metadata.csv
    ├── data/
    │   └── README.md
    ├── results/
    │   ├── tables/
    │   └── figures/
    └── docs/

Raw sequencing reads are not stored in this repository. Raw reads are/will be
available through the NCBI Sequence Read Archive under BioProject [ACCESSION].

## Requirements

Sequence processing was performed using:

- R [version]
- DADA2 [version]
- Cutadapt [version]
- [phyloseq, vegan, tidyverse, etc.]

See [session-info file/environment file] for complete package versions.

## Reproducibility

Scripts are numbered according to the order in which they should be run.

Raw FASTQ files should be placed in [directory] using the naming convention:

    SAMPLE_R1.fastq.gz
    SAMPLE_R2.fastq.gz

Run:

    scripts/01_cutadapt.sh

followed by:

    scripts/02_dada2.R

[Explain downstream workflow.]

## Data availability

Raw sequencing data will be deposited in the NCBI Sequence Read Archive (SRA)
upon publication. Processed ASV tables, taxonomy assignments, and sample
metadata are provided in this repository where appropriate.

## Citation

If you use data or code from this repository, please cite:

[Paper citation once available]

and/or:

[Zenodo DOI once repository is archived]

## Contact

Jacqueline Comstock
Marine Biological Laboratory
jcomstock@mbl.edu
