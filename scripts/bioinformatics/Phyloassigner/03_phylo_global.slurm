#!/bin/bash
#SBATCH --job-name=M197_global
#SBATCH --account=pi-azworden
#SBATCH --partition=caslake
#SBATCH --output=M197_global_%j.out
#SBATCH --error=M197_global_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

module load python/miniforge-25.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloassigner

cd /project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner

python /project/azworden/users/jcomstock/M197_V1V2/PhyloAssigner_local/pythonassigner_v0.9.py \
    -o 01_global \
    -r /project/azworden/databases/PhyloAssigner/databases/Global_16S_refDB/ref.aln \
    -t /project/azworden/databases/PhyloAssigner/databases/Global_16S_refDB/ref_tree.txt \
    -q /project/azworden/users/jcomstock/M197_V1V2/M197_V1V2_ASVs.fasta \
    -m /project/azworden/databases/PhyloAssigner/databases/Global_16S_refDB/edge.mapping \
    -T 8 \
    -p pplacer
