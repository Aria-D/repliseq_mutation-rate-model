# repliseq_mutation-rate-model

# HMEC Repli-seq Pipeline

This repository contains a Snakemake workflow to process Repli-seq data for HMEC cells, from raw FASTQs to per-base replication timing (RT) stored in a dense HDF5 format with G1, Early, and Late S-phase fractions.

---

## Project Structure

```markdown
project_root/
├── Snakefile
├── merge_repliseq_to_hdf5.py
├── config.yaml
├── samples.tsv
└── environment.yml
```
---

## Installation

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate repliseq_env
```

Ensure proper modules are loaded (if running on Slurm cluster):

```bash
module load bwa samtools bedtools
```

Configuration:

Edit the `config.yaml` file with the appropriate paths:

```yaml
samples: "samples.tsv"
genome: "genome.fa"
chrom_sizes: "genome.chrom.sizes"
outdir: "results"
hdf5: "HMEC_repliseq.h5"
```
- `samples.tsv` contains your sample metadata and FASTQ paths.

- `genome.fa` is the reference genome FASTA.

Make sure your genome file has an index built for bwa mem. If it does not, in the same directory as the reference genome run:

```bash
bwa index {/full/path/to/your_genome_file.fa}
```

Make sure you include the full path (`realpath`) to your directory that has your genome file!
