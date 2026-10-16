cd /project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner

awk -F'\t' '
NR > 1 && $4 ~ /^str(;|$)/ {print $3}
' 02_cyanoplastid/best_placements.tsv \
| sort -u > 02_cyanoplastid/stramenopile_ASV_ids.txt

module load python/miniforge-25.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloassigner

python - <<'PY'
from Bio import SeqIO

base = "/project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner"

with open(f"{base}/02_cyanoplastid/stramenopile_ASV_ids.txt") as f:
    ids = {line.strip() for line in f if line.strip()}

input_fasta = f"{base}/01_global/cyanoplastid_ASVs.fasta"
output_fasta = f"{base}/02_cyanoplastid/stramenopile_ASVs.fasta"

records = [
    record for record in SeqIO.parse(input_fasta, "fasta")
    if record.id in ids
]

SeqIO.write(records, output_fasta, "fasta")

print(f"Expected ASVs: {len(ids)}")
print(f"Written ASVs: {len(records)}")

if len(records) != len(ids):
    raise ValueError("ASV counts do not match!")
PY
