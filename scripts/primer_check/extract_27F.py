import csv

FASTA = "Prochlorococcus_full16S.fasta"
OUTPUT = "27F_sites.tsv"

primer = "AGAGTTTGATCMTGGCTCAG"

# Allowed bases for each primer position
allowed = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "M": {"A", "C"}
}

def read_fasta(filename):
    records = []
    header = None
    seq = []

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq).upper()))

                header = line[1:]
                seq = []

            else:
                seq.append(line)

        if header is not None:
            records.append((header, "".join(seq).upper()))

    return records


def mismatches(seq, primer):
    mismatch_positions = []

    for i, (base, p) in enumerate(zip(seq, primer), start=1):
        if base not in allowed[p]:
            mismatch_positions.append(i)

    return mismatch_positions


records = read_fasta(FASTA)

results = []

for header, sequence in records:

    best_site = None
    best_mismatches = None
    best_start = None

    for i in range(len(sequence) - len(primer) + 1):

        candidate = sequence[i:i + len(primer)]
        mm = mismatches(candidate, primer)

        if best_mismatches is None or len(mm) < len(best_mismatches):
            best_site = candidate
            best_mismatches = mm
            best_start = i + 1

    results.append({
        "header": header,
        "start": best_start,
        "site": best_site,
        "n_mismatches": len(best_mismatches),
        "mismatch_positions": ",".join(map(str, best_mismatches))
    })


with open(OUTPUT, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "header",
            "start",
            "site",
            "n_mismatches",
            "mismatch_positions"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(results)


print(f"Processed {len(results)} sequences")
print(f"Wrote results to {OUTPUT}")
