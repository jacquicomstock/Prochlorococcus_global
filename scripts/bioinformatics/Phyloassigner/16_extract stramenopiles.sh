python - <<'PY'
from Bio import SeqIO
import csv
from pathlib import Path

base = Path("/project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner")

placements = base / "04_stramenopiles/best_placements.tsv"
input_fasta = base / "02_cyanoplastid/stramenopile_ASVs.fasta"

groups = {
    "pel": set(),
    "dic": set()
}

with open(placements) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        taxon = row["taxon"].strip().rstrip(";")
        group = taxon.split(";")[0]
        if group in groups:
            groups[group].add(row["ASV_id"])

records = list(SeqIO.parse(input_fasta, "fasta"))

for group, ids in groups.items():
    output = base / f"04_stramenopiles/{group}_ASVs.fasta"
    selected = [r for r in records if r.id in ids]

    SeqIO.write(selected, output, "fasta")

    print(f"{group}: expected {len(ids)}, wrote {len(selected)}")

    if len(selected) != len(ids):
        raise ValueError(f"ASV mismatch for {group}")
PY
