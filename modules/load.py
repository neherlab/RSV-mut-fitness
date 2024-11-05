"""File to load and read counts for specific mutation types"""

import pandas as pd


def load_synonymous_muts(
    df_path,
    include_noncoding_muts=False,
    include_stop_tolerant_orfs=False,
    remove_orf9b=False,
):
    # Read in curated dataset
    df = pd.read_csv(df_path)

    # Create masks for synonymous mutations, non-coding mutations, stop-codon tolerant ORFs, and ORF9b
    mask_synonymous = (df["synonymous"] == True).values
    mask_noncoding = (df["noncoding"] == True).values
    stop_tolerant_orfs = ["ORF6", "ORF7a", "ORF7b", "ORF8", "ORF10"]
    pattern = "|".join(stop_tolerant_orfs)
    mask_tolerant = (df["gene"].str.contains(pattern)).values
    mask_orf9b = (
        (df["gene"] == "N;ORF9b")
        & (
            df["clade_founder_aa"].apply(lambda x: x[0])
            == df["mutant_aa"].apply(lambda x: x[0])
        )
    ).values

    # Apply masks
    if include_noncoding_muts:
        mask_synonymous = mask_synonymous + mask_noncoding
    if include_stop_tolerant_orfs:
        mask_synonymous = mask_synonymous + mask_tolerant
    if remove_orf9b:
        mask_synonymous = mask_synonymous + mask_orf9b

    return df[mask_synonymous]


def load_nonsynonymous_muts(df_path):
    # Read in curated dataset
    df = pd.read_csv(df_path)

    # Create masks for synonymous mutations, non-coding mutations, stop-codon tolerant ORFs, and ORF9b
    mask_nonsynonymous = (df.mut_class == "nonsynonymous").values

    return df[mask_nonsynonymous]


def load_nonexcluded_muts(df_path):
    # Read in curated dataset
    df = pd.read_csv(df_path)

    return df
