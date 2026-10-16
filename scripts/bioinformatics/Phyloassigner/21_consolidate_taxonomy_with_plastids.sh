cat > consolidate_taxonomy.py <<'PY'
import csv
from pathlib import Path
from collections import Counter

base = Path.cwd()

runs = {
    "global": "01_global",
    "cyanoplastid": "02_cyanoplastid",
    "cyanobacteria": "03_cyanobacteria",
    "stramenopiles": "04_stramenopiles",
    "pelagophyceae": "05_pelagophyceae",
    "dictyochophyceae": "06_dictyochophyceae",
    "viridiplantae": "07_viridiplantae",
}

# Ordered from broadest to most specialized.
priority = [
    "global",
    "cyanoplastid",
    "cyanobacteria",
    "stramenopiles",
    "viridiplantae",
    "pelagophyceae",
    "dictyochophyceae",
]

def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def clean_taxon(value):
    if not value or value.strip() in ("NA", ""):
        return ""
    return value.strip().rstrip(";")

data = {}
counts = {}

for name, folder in runs.items():
    best_file = base / folder / "best_placements.tsv"
    lca_file = base / folder / "LCA_placements.tsv"

    if not best_file.exists() or not lca_file.exists():
        raise FileNotFoundError(
            f"Missing classification files for {name}: {folder}"
        )

    best = read_tsv(best_file)
    lca = read_tsv(lca_file)

    best_ids = [r["ASV_id"] for r in best]
    lca_ids = [r["ASV_id"] for r in lca]

    if len(best_ids) != len(set(best_ids)):
        raise ValueError(f"Duplicate best ASV IDs in {name}")
    if len(lca_ids) != len(set(lca_ids)):
        raise ValueError(f"Duplicate LCA ASV IDs in {name}")
    if set(best_ids) != set(lca_ids):
        raise ValueError(f"Best/LCA ASV mismatch in {name}")

    counts[name] = len(best)

    for row in best:
        asv = row["ASV_id"]
        entry = data.setdefault(asv, {})
        entry[f"{name}_best"] = clean_taxon(row["taxon"])
        entry[f"{name}_LWR"] = row["like_weight_ratio"]
        entry[f"{name}_best_edge"] = row["edge_num"]

    for row in lca:
        asv = row["ASV_id"]
        entry = data.setdefault(asv, {})
        entry[f"{name}_LCA"] = clean_taxon(row["taxon"])
        entry[f"{name}_LCA_edge"] = row["edge_num"]

# Check specialized runs belong to their expected parent.
parents = {
    "cyanoplastid": "global",
    "cyanobacteria": "cyanoplastid",
    "stramenopiles": "cyanoplastid",
    "pelagophyceae": "stramenopiles",
    "dictyochophyceae": "stramenopiles",
    "viridiplantae": "cyanoplastid",
}

for child, parent in parents.items():
    child_ids = {
        asv for asv, row in data.items()
        if f"{child}_best" in row
    }
    parent_ids = {
        asv for asv, row in data.items()
        if f"{parent}_best" in row
    }
    missing = child_ids - parent_ids
    if missing:
        raise ValueError(
            f"{child}: {len(missing)} ASVs absent from {parent}"
        )

# Preserve full classifications; never overwrite parent
# taxonomy with an unresolved or missing specialized label.
for asv, entry in data.items():
    attempted = [
        name for name in priority
        if f"{name}_best" in entry
    ]
    entry["deepest_database_run"] = attempted[-1]

    for method in ("best", "LCA"):
        valid = [
            name for name in attempted
            if entry.get(f"{name}_{method}", "")
        ]
        if valid:
            selected = valid[-1]
            entry[f"final_{method}_taxonomy"] = (
                entry[f"{selected}_{method}"]
            )
            entry[f"final_{method}_database"] = selected
        else:
            entry[f"final_{method}_taxonomy"] = ""
            entry[f"final_{method}_database"] = ""

fields = ["ASV_id"]
for name in priority:
    fields.extend([
        f"{name}_best",
        f"{name}_LCA",
        f"{name}_LWR",
        f"{name}_best_edge",
        f"{name}_LCA_edge",
    ])

fields.extend([
    "deepest_database_run",
    "final_best_taxonomy",
    "final_best_database",
    "final_LCA_taxonomy",
    "final_LCA_database",
])

output = base / "M197_V1V2_PhyloAssigner_master_taxonomy.csv"

with open(output, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for asv in sorted(data):
        writer.writerow({"ASV_id": asv, **data[asv]})

print("\nCLASSIFICATION INPUTS")
for name in priority:
    print(f"{name}: {counts[name]} ASVs")

print(f"\nTotal unique ASVs: {len(data)}")
print(f"Master table: {output}")

print("\nFINAL BEST-PLACEMENT DATABASE COUNTS")
summary = Counter(
    r["final_best_database"] for r in data.values()
)
for name, count in summary.most_common():
    print(f"{name or 'unassigned'}: {count}")

print("\nFINAL LCA DATABASE COUNTS")
summary = Counter(
    r["final_LCA_database"] for r in data.values()
)
for name, count in summary.most_common():
    print(f"{name or 'unassigned'}: {count}")
PY
