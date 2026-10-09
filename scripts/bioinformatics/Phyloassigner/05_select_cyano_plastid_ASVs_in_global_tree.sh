OLD=/project/azworden/sequencing_archive/AVITI_AW-9031521_SIDE_A/M197/07_PhyloAssigner

grep "^>" $OLD/01_global/cyanoplastid_ASVs.fasta \
    | sed 's/^>//' \
    | sort -u \
    > /tmp/old_cyanoplastid_ids.txt

    awk -F'\t' '
NR==FNR {x[$1]=1; next}
NR>1 && ($3 in x) {print $1}
' \
/tmp/old_cyanoplastid_ids.txt \
$OLD/01_global/best_placements.tsv \
| sort -u \
> /tmp/cyanoplastid_edges.txt

awk -F'\t' '
NR==FNR {edges[$1]=1; next}
NR>1 && ($1 in edges) {print $3}
' \
/tmp/cyanoplastid_edges.txt \
07_PhyloAssigner/01_global/best_placements.tsv \
| sort -u \
> /tmp/new_cyanoplastid_ids.txt

python - <<'PY'
from Bio import SeqIO

ids = set(
    line.strip()
    for line in open("/tmp/new_cyanoplastid_ids.txt")
)

input_fasta = (
    "/project/azworden/users/jcomstock/"
    "M197_V1V2/M197_V1V2_ASVs.fasta"
)

output_fasta = (
    "/project/azworden/users/jcomstock/M197_V1V2/"
    "07_PhyloAssigner/01_global/cyanoplastid_ASVs.fasta"
)

records = [
    r for r in SeqIO.parse(input_fasta, "fasta")
    if r.id in ids
]

SeqIO.write(
    records,
    output_fasta,
    "fasta"
)

print(f"Wrote {len(records)} sequences")
print(f"Expected {len(ids)} sequences")
PY
