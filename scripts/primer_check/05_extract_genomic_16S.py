import re
import os

INPUT_DIR = "strain_genbank"
OUTPUT = "cultured_strain_16S_v2.fasta"


def reverse_complement(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def get_genome_sequence(lines):
    in_origin = False
    seq = []

    for line in lines:
        if line.startswith("ORIGIN"):
            in_origin = True
            continue

        if in_origin:
            if line.startswith("//"):
                break

            seq.append("".join(re.findall(r"[acgtACGT]", line)))

    return "".join(seq).upper()


records = []

for filename in sorted(os.listdir(INPUT_DIR)):

    if not filename.endswith(".gb"):
        continue

    strain = filename.replace(".gb", "")
    path = os.path.join(INPUT_DIR, filename)

    with open(path) as f:
        lines = f.readlines()

    genome = get_genome_sequence(lines)

    copy_number = 0

    for i, line in enumerate(lines):

        if not line.startswith("     rRNA"):
            continue

        # Find the end of THIS feature.
        j = i + 1

        while j < len(lines):

            # A new GenBank feature begins when columns 6-20
            # contain a feature name.
            if (
                lines[j].startswith("     ")
                and len(lines[j]) > 20
                and lines[j][5:21].strip()
                and not lines[j][5:21].lstrip().startswith("/")
            ):
                break

            j += 1

        feature_block = "".join(lines[i:j])

        # Only accept this exact rRNA feature if its own
        # annotation identifies it as 16S.
        if '/product="16S ribosomal RNA' not in feature_block:
            continue

        location = line[21:].strip()

        copy_number += 1

        is_complement = location.startswith("complement")

        nums = re.findall(r"\d+", location)

        start = int(nums[0])
        end = int(nums[-1])

        seq = genome[start - 1:end]

        if is_complement:
            seq = reverse_complement(seq)

        records.append(
            (strain, copy_number, location, seq)
        )


with open(OUTPUT, "w") as out:

    for strain, copy_number, location, seq in records:

        header = (
            f">{strain}_16S_copy{copy_number} "
            f"location={location}"
        )

        out.write(header + "\n")

        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + "\n")


print("Extracted", len(records), "16S sequences")
print("Wrote", OUTPUT)

for strain, copy_number, location, seq in records:

    print(
        strain,
        "copy", copy_number,
        "length", len(seq),
        "location", location
    )
