"""
Script for computing fitness effects of nucleotide mutations

Input:
- `muts_clade_counts.csv` --> table with clade-wise predicted counts
- `clade_clusteing.json` --> dicitionary for grouping clades

Output:
- `{group}_nmut_fitness` --> table with fitness effects for each group

"""

# Required packages
import pandas as pd
import sys
import os
import json

# Adding module folder to system paths
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

# Import module for probabilistic fitness estimates
from modules import probfit

# Read-in clades grouping
with open(sys.argv[1]) as json_clade:
    clade_cluster = json.load(json_clade)

# Import table with clade-wise counts
muts_by_clade = pd.read_csv("../results/mut_counts_by_clade.csv", low_memory=False)

# Columns not be aggregated
group_cols = [
    "nt_mutation",
    "gene",
    "codon_site",
    "aa_mutation",
    "synonymous",
    "noncoding",
]

# Compute fitness for each group and save to file
for key, val in clade_cluster.items():
    muts_by_clade_cluster = (
        muts_by_clade.query("clade.isin(@val)")
        .groupby(group_cols, as_index=False)
        .aggregate(
            expected_count=pd.NamedAgg("expected_count", "sum"),
            predicted_count=pd.NamedAgg("predicted_count", "sum"),
            actual_count=pd.NamedAgg("actual_count", "sum"),
            tau_squared=pd.NamedAgg("tau_squared", "sum"),
        )
    )
    muts_by_clade_cluster.insert(
        0, "nt_site", muts_by_clade_cluster["nt_mutation"].apply(lambda x: int(x[1:-1]))
    )
    muts_by_clade_cluster.insert(0, "cluster", key)
    muts_by_clade_cluster = muts_by_clade_cluster.sort_values("nt_site").reset_index(
        drop=True
    )
    probfit.add_probabilistic_estimates(muts_by_clade_cluster, N_f=300)
    muts_by_clade_cluster.to_csv(
        f"../results/ntmut_fitness/{key}_ntmut_fitness.csv", index=False
    )
