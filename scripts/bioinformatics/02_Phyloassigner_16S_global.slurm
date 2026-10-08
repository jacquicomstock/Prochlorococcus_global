#!/bin/bash
#SBATCH --job-name=M197_phylo_global
#SBATCH --account=pi-azworden
#SBATCH --partition=caslake
#SBATCH --output=M197_phylo_global_%j.out
#SBATCH --error=M197_phylo_global_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /project/azworden/databases/PhyloAssigner

module load python/miniforge-25.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloassigner

echo "Job started: $(date)"
echo "Running on: $(hostname)"

python pythonassigner_v0.9.py \
    --out_dir /project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner/01_global \
    --ref_align /project/azworden/databases/PhyloAssigner/databases/Global_16S_refDB/ref.aln \
    --ref_tree /project/azworden/databases/PhyloAssigner/databases/Global_16S_refDB/ref_tree.txt \
    --query_seqs /project/azworden/users/jcomstock/M197_V1V2/M197_V1V2_ASVs.fasta \
    --mapping /project/azworden/databases/PhyloAssigner/databases/Global_16S_refDB/edge.mapping \
    --threads 8 \
    --placer pplacer
