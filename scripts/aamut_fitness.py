"""
Script for estimating amino acid fitness for all groups of clades
"""

# Import packages
import pandas as pd
import sys
import os

# Adding module folder to system paths
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

from modules import aamutfit

# How to handle nucleotide sites that overlap among sites: exclude them from
# amino-acid fitness estimates or retain those sites in estimates. If retained, estimates
# can be confounded by not knowing which gene the selection is acting on.
gene_overlaps = {
    "exclude": [["ORF1a", "ORF1ab"], ["ORF7a", "ORF7b"]],
    "retain": [["N", "ORF9b"]],
}

# Orf1ab to Nsp numbering (amino-acid start in Orf1ab) from
# https://github.com/theosanderson/Codon2Nucleotide/blob/main/src/App.js
orf1ab_to_nsps = {
    "nsp1": [1, 180],
    "nsp2": [181, 818],
    "nsp3": [819, 2763],
    "nsp4": [2764, 3263],
    "nsp5 (Mpro)": [3264, 3569],
    "nsp6": [3570, 3859],
    "nsp7": [3860, 3942],
    "nsp8": [3943, 4140],
    "nsp9": [4141, 4253],
    "nsp10": [4254, 4392],
    "nsp12 (RdRp)": [4393, 5324],
    "nsp13": [5325, 5925],
    "nsp14": [5926, 6452],
    "nsp15": [6453, 6798],
    "nsp16": [6799, 7096],
}

explode_cols = [
    "gene",
    "clade_founder_aa",
    "mutant_aa",
    "codon_site",
    "aa_mutation",
]

# Fitness pseudocount
fitness_pseudocount = 0.5

# Directory containing nucleotide fitness files
nmut_fit_directory = "../results/ntmut_fitness/"

for name in os.listdir(nmut_fit_directory):
    # Open file
    with open(os.path.join(nmut_fit_directory, name)) as f:
        # Read content of file
        nmut_fit = pd.read_csv(f)
        cluster = nmut_fit.loc[0, "cluster"]
        # Get only coding mutations
        nmut_fit_coding = aamutfit.get_coding(nmut_fit, gene_overlaps, explode_cols)
        # Aggregate counts for amino acid mutations
        aa_counts = aamutfit.aggregate_counts(nmut_fit_coding, explode_cols)
        # Adding naive fitness estimates to `aa_counts`
        aamutfit.naive_fitness(aa_counts, fitness_pseudocount=fitness_pseudocount)
        # Dataframe with refined fitness estimates
        aa_fit = aamutfit.aa_fitness(nmut_fit_coding, explode_cols)
        # Define dataframe for ORF1ab to nsps mapping
        orf1ab_to_nsps_df = aamutfit.map_orf1ab_to_nsps(orf1ab_to_nsps)
        # Adding rows for nsps
        aa_counts = aamutfit.add_nsps(aa_counts, orf1ab_to_nsps_df)
        aa_fit = aamutfit.add_nsps(aa_fit, orf1ab_to_nsps_df)
        # Merging counts and fitness dataframes
        aamut_fitness = aamutfit.merge_aa_df(aa_fit, aa_counts, explode_cols)
        # Save amino acid mutations fitness to file
        aamut_fitness.to_csv(
            f"../results/aamut_fitness/{cluster}_aamut_fitness.csv", index=False
        )
