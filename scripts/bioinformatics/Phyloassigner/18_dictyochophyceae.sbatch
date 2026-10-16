cat > M197_dictyochophyceae.sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=M197_dictyochophyceae
#SBATCH --account=pi-azworden
#SBATCH --partition=caslake
#SBATCH --output=M197_dictyochophyceae_%j.out
#SBATCH --error=M197_dictyochophyceae_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

module load python/miniforge-25.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloassigner

cd /project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner

python /project/azworden/users/jcomstock/M197_V1V2/PhyloAssigner_local/pythonassigner_v0.9.py \
    -o 06_dictyochophyceae \
    -r /project/azworden/databases/PhyloAssigner/databases/Dictyochophyceae_refDB/ref.aln \
    -t /project/azworden/databases/PhyloAssigner/databases/Dictyochophyceae_refDB/ref_tree.txt \
    -q 04_stramenopiles/dic_ASVs.fasta \
    -m updated_mappings/dictyo.edge.mapping.updated10.26 \
    -T 8 \
    -p pplacer
EOF
