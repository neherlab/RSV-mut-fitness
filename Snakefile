""" Snakemake file for running the pipeline """

import yaml

configfile: "config.yaml"

rule all:
    input:
        expand('results/aamut_fitness/{cluster}_aamut_fitness.csv', cluster=config['clade_cluster'].keys()),

rule get_counts_table:
    params:
        url_counts=config["url_counts"],
    output:
        csv=temp('results/expected_vs_actual_counts.csv')
    shell:
        """
            curl -k {params.url_counts} > {output.csv}
        """
rule get_clade_founder:
    params:
        url_founder=config["url_founder"],
    output:
        csv=temp('results/clade_founder.csv')
    shell:
        """
            curl {params.url_founder} > {output.csv}
        """

rule annotate_counts:
    input:
        rna_struct="data/lan_2022/41467_2022_28603_MOESM11_ESM.txt",
        counts=rules.get_counts_table.output.csv,
        clade_founder=rules.get_clade_founder.output.csv,
    output:
        csv=temp("results/mut_counts_by_clade.csv"),
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

rule predicted_counts:
    input:
        counts_df=rules.annotate_counts.output.csv,
        pre_omicron=rules.curated_counts.output.pre_omicron,
        omicron=rules.curated_counts.output.omicron,
    output:
        pred_count_csv=temp("results/pred_mut_counts_by_clade.csv"),
    notebook:
        "notebook/predicted_counts_by_clade.py.ipynb"

rule counts_cluster:
    params:
        cluster=lambda wc: wc.cluster,
        clades=lambda wc: config['clade_cluster'][wc.cluster],
    input:
        counts_df=rules.predicted_counts.output.pred_count_csv,
    output:
        cluster_counts=temp('results/ntmut_fitness/{cluster}_ntmut_counts.csv'),
    notebook:
        "notebook/ntmut_counts_cluster.py.ipynb"

rule ntmut_fitness:
    input:
        cluster_counts='results/ntmut_fitness/{cluster}_ntmut_counts.csv' 
    output:
        ntfit_csv='results/ntmut_fitness/{cluster}_ntmut_fitness.csv',
    notebook:
        'notebook/ntmut_fitness.py.ipynb'

rule aamut_fitness:
    params:
        orf_to_nsps=config['orf1ab_to_nsps'],
        gene_ov=config['gene_overlaps'],
        fit_pseudo=config['fitness_pseudocount'],
    input:
        ntfit_csv='results/ntmut_fitness/{cluster}_ntmut_fitness.csv',
    output:
        aafit_csv='results/aamut_fitness/{cluster}_aamut_fitness.csv',
    notebook:
        'notebook/aamut_fitness.py.ipynb'

localrules:
    get_counts_table,
    get_clade_founder,
