##########################################
# Snakefile: Repli-seq processing pipeline
##########################################

import os
import pandas as pd
import h5py
from functools import reduce

SAMPLES = config["samples"]
GENOME = config["genome"]
OUTDIR = config["outdir"]
HDF5_NAME = config["hdf5"]

# make sure output directory exists
os.makedirs(OUTDIR, exist_ok=True)

rule all:
    input:
        os.path.join(OUTDIR, HDF5_NAME)

# load sample sheet
samples = pd.read_csv(SAMPLES, sep="\t")

import pandas as pd

samples = pd.read_csv("samples.tsv", sep="\t")
samples["Sample"] = samples["Sample"].str.strip()
samples["FQ1"] = samples["FQ1"].str.strip()
samples["FQ2"] = samples["FQ2"].str.strip()

def get_fastqs(wildcards):
    row = samples.loc[samples["Sample"] == wildcards.sample]
    print(row)
    if row.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in samples.tsv")
    print([row["FQ1"].values[0], row["FQ2"].values[0]])
    return [row["FQ1"].values[0], row["FQ2"].values[0]]


rule trim:
    input: get_fastqs
    output:
        R1=os.path.join(OUTDIR, "{sample}_R1.trim.fastq.gz"),
        R2=os.path.join(OUTDIR, "{sample}_R2.trim.fastq.gz")
    shell:
        """
        cutadapt -q 0 -O 1 -m 0 -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
        -o {output.R1} -p {output.R2} {input[0]} {input[1]}
        """

rule align:
    input:
        r1=os.path.join(OUTDIR, "{sample}_R1.trim.fastq.gz"),
        r2=os.path.join(OUTDIR, "{sample}_R2.trim.fastq.gz")
    output:
        bam=os.path.join(OUTDIR, "{sample}.bam")
    threads: 8
    shell:
        """
        bwa mem -t {threads} {GENOME} {input.r1} {input.r2} | \
        samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule filter_bam:
    input:
        bam=os.path.join(OUTDIR, "{sample}.bam")
    output:
        bam=os.path.join(OUTDIR, "{sample}.filtered.bam"),
        bai=os.path.join(OUTDIR, "{sample}.filtered.bam.bai")
    shell:
        """
        samtools view -b -q 30 -F 4 {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule coverage:
    input:
        bam=os.path.join(OUTDIR, "{sample}.filtered.bam")
    output:
        cov=os.path.join(OUTDIR, "{sample}.perbase.cov")
    shell:
        """
        bedtools genomecov -d -ibam {input.bam} > {output.cov}
        """

rule merge_to_hdf5:
    input:
        covs=expand(os.path.join(OUTDIR, "{sample}.perbase.cov"), sample=samples.Sample)
    output:
        h5=os.path.join(OUTDIR, HDF5_NAME)
    shell:
        """
        python merge_repliseq_to_hdf5.py {config[samples]} {config[genome]} {config[outdir]} {config[hdf5]}
        """

