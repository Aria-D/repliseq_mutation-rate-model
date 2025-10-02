# repliseq_mutation-rate-model

# HMEC Repli-seq Pipeline

This repository contains a Snakemake workflow to process Repli-seq data for HMEC cells, from raw FASTQs to per-base replication timing (RT) stored in a dense HDF5 format with G1, Early, and Late S-phase fractions.

---

## Project Structure

project_root/
├── Snakefile
├── merge_repliseq_to_hdf5.py
├── config.yaml
├── samples.tsv
└── environment.yml
---

## Installation

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate repliseq_env

Ensure proper modules are loaded (if running on Slurm cluster):

```bash
module load bwa samtools bedtools

