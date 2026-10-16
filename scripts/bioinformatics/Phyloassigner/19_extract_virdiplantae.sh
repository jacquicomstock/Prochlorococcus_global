cd /project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner

module load python/miniforge-25.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phyloassigner

python - <<'PY'
from Bio import SeqIO
import csv
from pathlib import Path

base = Path("/project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner")

placements = base / "02_cyanoplastid/best_placements.tsv"
input_fasta = base / "01_global/cyanoplastid_ASVs.fasta"
output_fasta = base / "02_cyanoplastid/vir_ASVs.fasta"

ids = set()

with open(placements) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        taxon = row["taxon"].strip().rstrip(";")
        if taxon == "vir" or taxon.startswith("vir;"):
            ids.add(row["ASV_id"])

records = [
    r for r in SeqIO.parse(input_fasta, "fasta")
    if r.id in ids
]

SeqIO.write(records, output_fasta, "fasta")

print(f"Expected ASVs: {len(ids)}")
print(f"Written ASVs: {len(records)}")

if len(records) != len(ids):
    raise ValueError("ASV count mismatch!")
PY
