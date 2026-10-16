cat > M197_stramenopiles.sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=M197_stramenopiles
#SBATCH --account=pi-azworden
#SBATCH --partition=caslake
#SBATCH --output=M197_stramenopiles_%j.out
#SBATCH --error=M197_stramenopiles_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

module load python/miniforge-25.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloassigner

cd /project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner

python /project/azworden/users/jcomstock/M197_V1V2/PhyloAssigner_local/pythonassigner_v0.9.py \
    -o 04_stramenopiles \
    -r /project/azworden/databases/PhyloAssigner/databases/Stramenopiles_refDB/ref.aln \
    -t /project/azworden/databases/PhyloAssigner/databases/Stramenopiles_refDB/ref_tree.txt \
    -q 02_cyanoplastid/stramenopile_ASVs.fasta \
    -m /project/azworden/databases/PhyloAssigner/databases/Stramenopiles_refDB/edge.mapping \
    -T 8 \
    -p pplacer
EOF
