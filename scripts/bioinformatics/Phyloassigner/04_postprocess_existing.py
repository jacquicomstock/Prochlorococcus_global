import os
import pandas as pd

from modules.best import getBEST
from modules.LCA import getLCA
from modules.helper_functions import merge_mapping

out_dir = "/project/azworden/users/jcomstock/M197_V1V2/07_PhyloAssigner/01_global"

placements = os.path.join(
    out_dir,
    "placements.jplace"
)

mapping_file = (
    "/project/azworden/databases/PhyloAssigner/"
    "databases/Global_16S_refDB/edge.mapping"
)

# Generate best placements
best = getBEST(
    placements,
    "pplacer"
)

mapping = pd.read_csv(
    mapping_file,
    sep="\t",
    header=None
).set_axis(
    ["edge_num", "taxon"],
    axis=1
)

best = merge_mapping(
    best,
    mapping
)

best.to_csv(
    os.path.join(out_dir, "best_placements.tsv"),
    sep="\t",
    index=False
)

# Generate LCA placements
getLCA(
    out_dir,
    placements,
    "8"
)
