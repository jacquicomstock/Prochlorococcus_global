# the purpose of this script is to check the surrounding 27F site to see whether 
# we can trust the variations in 27F primer site sequences that we observed.

FASTA = "Prochlorococcus_full16S_aligned.fasta"

targets = [
    "DQ070784.1",
    "KJ589672.1",
    "KU498065.1",
    "KX179943.1",
    "KX179953.1",
    "KX987406.1",
    "KX987524.1",
    "KY561501.1",
    "KY561511.1",
    "MF488575.1",
    "MF540170.1"
]

REFERENCE_SITE = "AGAGTTTGATCCTGGCTCAG"


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

# Find a reference containing the canonical site
for header, seq in records:
    ungapped = seq.replace("-", "")
    pos = ungapped.find(REFERENCE_SITE)

    if pos != -1:
        ref = seq

        # Determine alignment columns for start/end of site
        ungapped_pos = 0
        cols = []

        for col, base in enumerate(ref):
            if base != "-":
                if pos <= ungapped_pos < pos + len(REFERENCE_SITE):
                    cols.append(col)
                ungapped_pos += 1

        start = cols[0]
        end = cols[-1]
        break


print("27F alignment columns:", start + 1, "-", end + 1)
print()

for target in targets:
    for header, seq in records:
        if header.startswith(target):

            left = max(0, start - 30)
            right = min(len(seq), end + 31)

            context = seq[left:right]

            print(">", header)
            print(context)
            print(" " * (start-left) + "^" * (end-start+1))
            print()
