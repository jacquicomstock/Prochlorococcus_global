import csv

FASTA = "Prochlorococcus_full16S_aligned.fasta"
OUTPUT = "27F_sites_aligned.tsv"

# A known perfect Prochlorococcus 27F-site variant
REFERENCE_SITE = "AGAGTTTGATCCTGGCTCAG"

PRIMER = "AGAGTTTGATCMTGGCTCAG"

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


records = read_fasta(FASTA)

# Find a sequence containing the known 27F site after removing gaps
reference_header = None
reference_aligned = None

for header, aligned_seq in records:
    ungapped = aligned_seq.replace("-", "")
    if REFERENCE_SITE in ungapped:
        reference_header = header
        reference_aligned = aligned_seq
        break

if reference_aligned is None:
    raise ValueError("Could not find a sequence containing the reference 27F site.")


# Find alignment columns corresponding to the 21 nt reference site
ungapped_pos = 0
site_start_ungapped = reference_aligned.replace("-", "").find(REFERENCE_SITE)

site_columns = []

for col, base in enumerate(reference_aligned):
    if base != "-":
        if site_start_ungapped <= ungapped_pos < site_start_ungapped + len(REFERENCE_SITE):
            site_columns.append(col)
        ungapped_pos += 1


print("Reference sequence:", reference_header)
print("27F alignment columns:", site_columns[0] + 1, "to", site_columns[-1] + 1)


results = []

for header, aligned_seq in records:

    site = "".join(aligned_seq[col] for col in site_columns)

    # If sequence does not cover this region, report NA
    if set(site) == {"-"}:
        status = "no_coverage"
        n_mismatches = "NA"
        mismatch_positions = "NA"

    elif "-" in site:
        status = "partial_or_gap"
        n_mismatches = "NA"
        mismatch_positions = "NA"

    else:
        mismatch_positions_list = []

        for i, (base, p) in enumerate(zip(site, PRIMER), start=1):
            if base not in allowed[p]:
                mismatch_positions_list.append(i)

        status = "complete"
        n_mismatches = len(mismatch_positions_list)
        mismatch_positions = ",".join(map(str, mismatch_positions_list))

    results.append({
        "header": header,
        "site": site,
        "status": status,
        "n_mismatches": n_mismatches,
        "mismatch_positions": mismatch_positions
    })


with open(OUTPUT, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "header",
            "site",
            "status",
            "n_mismatches",
            "mismatch_positions"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(results)


print("Processed", len(results), "sequences")
print("Wrote", OUTPUT)
