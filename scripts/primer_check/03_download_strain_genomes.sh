#!/bin/bash

mkdir -p strain_genomes_v3

while IFS= read -r strain <&3
do
    echo "Searching for $strain..."

    esearch -db nucleotide \
        -query "\"$strain\"[strain] AND Prochlorococcus[Organism] AND \"complete genome\"[Title]" \
        | efetch -format fasta \
        > "strain_genomes_v3/${strain}.fasta"

    n=$(grep -c '^>' "strain_genomes_v3/${strain}.fasta")

    echo "  Found $n record(s)"

done 3< Prochlorococcus_strains.txt
