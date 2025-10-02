#!/bin/bash
#SBATCH --job-name=repliseq_pipeline
#SBATCH --output=logs/repliseq_%j.out
#SBATCH --error=logs/repliseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=scu-cpu     # change if needed

# Load modules (if needed)
module load bwa
module load samtools
module load bedtools

# Activate conda environment with snakemake + python packages
conda activate repliseq_env

# Make logs folder if not exists
mkdir -p logs

echo "[INFO] Starting Repli-seq Snakemake pipeline on $(date)"

# Run Snakemake
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores $SLURM_CPUS_PER_TASK \
    --use-conda \
    --rerun-incomplete \
    --printshellcmds \
    --latency-wait 60 \
    --directory .

echo "[INFO] Finished pipeline on $(date)"
