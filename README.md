# mouseCRC

This repository contains the anlysis code for the paper "Spatially defined multicellular functional units in colorectal cancer revealed from single cell and spatial transcriptomics" https://doi.org/10.1101/2022.10.02.508492 .

## How to use this repository

The data processing is set up as a snakemake (https://snakemake.readthedocs.io/en/stable/) workflow. It does not need to be installed as it only describes the workflow to be run by snakemake. It should run on any current Linux (tested on RHEL 7.9, with snakemake 7.19).

After setting up snakemake, the workflow can be started from the main folder using the command `snakemake`. For more information about running workflows with snakemake, see the snakemake docs (https://snakemake.readthedocs.io/en/stable/). In particular, it is recommended to run it in an HPC cluster environment with cluster specific snakemake profile. Running on a standard desktop computer is expected to take days if not weeks.

The workflow will create all necessary environments, download publically available datasets, populate a new directory `notebooks` with all the executed jupyter notebooks for plot generation, and collect all plots in a new `results/figures` directory.


