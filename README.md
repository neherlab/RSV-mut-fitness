# Fitness effects of nucleotide and amino acid mutations updated with novel estimates of neutral rates.

## Overview
This repository computes an updated estimate of the fitness effects of mutations of the SARS-CoV-2 genome,
expanding and building on the approach presented in a [paper](https://academic.oup.com/ve/article/9/2/vead055/7265011)
by [Jesse Bloom](https://scholar.google.com/citations?user=S12x_eQAAAAJ&hl=en) and [Richard Neher](https://neherlab.org/).

The two main novelties consist in:
* Estimating the mutation rates with a *general linear model* that takes into account several factors such as:
  - A site's local nucleotide context.
  - Its structural pairing in RNA secondary structure.
  - The region of the genome it belongs to.
* Fitness effects of mutations are now computed within a Bayesian probabilistic framework that also provides uncertainties.

## References
- Details about the computational framework can be found  in the related [paper](https://github.com/matsengrp/SARS2-synonymous-mut-rate-tex).
- Data used as reference for RNA secondary structure are in [Lan et al](https://www.nature.com/articles/s41467-022-28603-2).
- Evidence about influence of secondary structure on mutation rates were first presented in a paper by [Hensel](https://www.biorxiv.org/content/10.1101/2024.02.27.581995v1.abstract).
- The original approach for estimating mutational fitness is presented in [Bloom & Neher](https://academic.oup.com/ve/article/9/2/vead055/7265011).

## Computational pipeline
It is possible to reproduce the fitness estimates by running the computational analysis defined in this GitHub repository.

Firstly, a customized [conda](https://docs.conda.io/) environment needs to be buily from the [environment.yml](environment.yml) file. To do so, you must install [conda](https://docs.conda.io/) and then run:

    conda env create -f environment.yml

This will create a [conda](https://docs.conda.io/) environment called `SARS2-refined-fitness`, that you need to activate:

    conda activate SARS2-refined-fitness

The pipeline is managed by [snakemake](https://snakemake.readthedocs.io/) through a [Snakefile](Snakefile), whose configuration is defined in [config.yaml](config.yaml). To run the pipeline use:

    snakemake -c <n_cpus>

where `n_cpus` is the number of CPU cores you want to use.

The pipeline mainly relies on Python Jupyter notebooks to run. These can be found in the [./notebook](notebook) folder.

### Input

The pipeline runs downstream from two files fetched directly from the [jbloomlab/SARS2-mut-fitness](https://github.com/jbloomlab/SARS2-mut-fitness/tree/main) GitHub repo:

- Table of clade founders nucleotides [~/results_gisaid_2024-04-24/clade_founder_nts/clade_founder_nts.csv](https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/results_gisaid_2024-04-24/clade_founder_nts/clade_founder_nts.csv).
- Clade-wise table of counts [~/results_gisaid_2024-04-24/expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv](https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/results_gisaid_2024-04-24/expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv).

The related links are defined in the [config.yaml](config.yaml) file and can be changed to any version of the reference UShER tree.

### Configuration
Ahead of the computation of mutational fitness effects, predicted and actual mutation counts can be aggregated by defining clusters of clades. This is defined by a dictionary `clade_clusters` in the [config.yaml](config.yaml) file, which can be customized.

### Output
Files produced by the pipeline are saved into the [./results](results) folder. These are subsequently divided in the following subfolders:
- [curated](results/curated/): the training datasets for the *General Linear Model* (GLM) producing the predicted counts. These are divided among pre and Omicron clades.
- [master_tables](results/master_tables/): the reference models inferred from the training datasets, i.e. a site's context and the associated mutation rate.
- [ntmut_fitness](results/ntmut_fitness/): files `{cluster}_ntmut_fitness.csv` for each cluster of clades with the nucleotide mutation fitness effects.
- [aamut_fitness](results/aamut_fitness/): files `{cluster}_aamut_fitness.csv` with fitness effects of the amino acid mutations for each cluster of clades.
