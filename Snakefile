"""
This pipeline fetches raw mutation counts from Jesse Bloom's GitHub repository
and produces refined fitness including estimates of uncertainty for each mutation
in different subsets of the total sequence availability.
"""

configfile: "config.yaml"

rule all:
    input:
        #expand('results/aamut_fitness/{cluster}_aamut_fitness.csv', cluster=config['clade_cluster'].keys()),
        'results/aamut_fitness/aamut_fitness_by_cluster.csv'

rule get_counts_table:
    message:
        "Downloading table with mutation counts from Jesse Bloom's GitHub repository"
    params:
        url_counts=config["url_counts"],
    output:
        csv=temp('results/expected_vs_actual_counts.csv')
    shell:
        """
            curl -k {params.url_counts} > {output.csv}
        """
rule get_clade_founder:
    message:
        "Downloading table with clade founder sequences from Jesse Bloom's GitHub repository"
    params:
        url_founder=config["url_founder"],
    output:
        csv=temp('results/clade_founder.csv')
    shell:
        """
            curl {params.url_founder} > {output.csv}
        """

rule annotate_counts:
    message:
        "Augment the counts table with RNA secondary structure pairing information, sequence context, and other features"
    input:
        rna_struct="data/lan_2022/41467_2022_28603_MOESM11_ESM.txt",
        counts=rules.get_counts_table.output.csv,
        clade_founder=rules.get_clade_founder.output.csv,
    output:
        csv=temp("results/mut_counts_by_clade.csv"),
    notebook:
        "notebook/counts_by_clade.py.ipynb"

rule curated_counts:
    message:
        "Create training dataset to infer the General Linear Model for mutations in Omicron and pre-Omicron sequences"
    input:
        mut_counts=rules.annotate_counts.output.csv,
    output:
        pre_omicron="results/curated/curated_mut_counts_pre_omicron.csv",
        omicron="results/curated/curated_mut_counts_omicron.csv"
    notebook:
        "notebook/curate_counts_pre_post_omicron.py.ipynb"

rule master_table:
    message:
        "Create tables with predicted mutation rates for each mutation in each of its contexts for pre-Omicron and Omicron sequences"
    input:
        pre_omicron_counts=rules.curated_counts.output.pre_omicron,
        omicron_counts=rules.curated_counts.output.omicron,
    output:
        pre_omicron_ms='results/master_tables/master_table_pre_omicron.csv',
        omicron_ms='results/master_tables/master_table_omicron.csv',
    notebook:
        "notebook/master_tables.py.ipynb"

rule predicted_counts:
    message:
        "Add predicted counts, based on the inferred mutation rate model, to the table with observed mutation counts."
    input:
        counts_df=rules.annotate_counts.output.csv,
        pre_omicron=rules.curated_counts.output.pre_omicron,
        omicron=rules.curated_counts.output.omicron,
    output:
        pred_count_csv=temp("results/pred_mut_counts_by_clade.csv"),
    notebook:
        "notebook/predicted_counts_by_clade.py.ipynb"

rule counts_cluster:
    message:
        """
        Create tables for each subset of sequences (clades, groups of clades, etc.) that contain the actual and prediced counts.
        These groups are defined in the config file as 'clade_clusters'.
        """
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
    message:
        "Calculate estimates of the fitness effects of each nucleotide mutation."
    input:
        cluster_counts='results/ntmut_fitness/{cluster}_ntmut_counts.csv'
    output:
        ntfit_csv='results/ntmut_fitness/{cluster}_ntmut_fitness.csv',
    notebook:
        'notebook/ntmut_fitness.py.ipynb'

rule aamut_fitness:
    message:
        "Calculate estimates of the fitness effects of each amino acid substitution."
    params:
        orf_to_nsps=config['orf1ab_to_nsps'],
        gene_ov=config['gene_overlaps'],
        genes=config['genes'],
        fit_pseudo=config['fitness_pseudocount'],
    input:
        ntfit_csv='results/ntmut_fitness/{cluster}_ntmut_fitness.csv',
    output:
        aafit_csv='results/aamut_fitness/{cluster}_aamut_fitness.csv',
    notebook:
        'notebook/aamut_fitness.py.ipynb'

rule concat_aamut:
    message:
        "Concatenating {{cluster_aamut}}_fitness.csv files"
    input:
        aafit_csv=expand('results/aamut_fitness/{cluster}_aamut_fitness.csv', cluster=config['clade_cluster'].keys())
    output:
        aafit_concat=temp('results/aamut_fitness/aamut_fitness_by_cluster.csv')
    shell:
        """
        {{ 
            head -n 1 {input.aafit_csv[0]};
            tail -n +2 -q {input.aafit_csv}
        }} > {output.aafit_concat}
        """

