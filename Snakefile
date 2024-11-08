""" Snakemake file for running the pipeline """

import yaml

configfile: "config.yaml"

rule all:
    input:
        'results/master_tables/master_table_pre_omicron.csv',
        'results/master_tables/master_table_omicron.csv',

rule annotate_counts:
    params:
        url_founder=config["url_founder"],
        url_counts=config["url_counts"],
    input:
        rna_struct="data/lan_2022/41467_2022_28603_MOESM11_ESM.txt"
    output:
        csv="results/mut_counts_by_clade.csv",
    notebook:
        "notebook/counts_by_clade.py.ipynb"

rule curated_counts:
    input:
        mut_counts=rules.annotate_counts.output.csv,
    output:
        pre_omicron="results/curated/curated_mut_counts_pre_omicron.csv",
        omicron="results/curated/curated_mut_counts_omicron.csv"
    notebook:
        "notebook/curate_counts_pre_post_omicron.py.ipynb"

rule master_table:
    input:
        pre_omicron_counts=rules.curated_counts.output.pre_omicron,
        omicron_counts=rules.curated_counts.output.omicron,
    output:
        pre_omicron_ms='results/master_tables/master_table_pre_omicron.csv',
        omicron_ms='results/master_tables/master_table_omicron.csv',
    notebook:
        "notebook/master_tables.py.ipynb"