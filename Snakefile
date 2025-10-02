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

rule trim:
    input:
        fq1=lambda wildcards: samples.loc[samples.Sample==wildcards.sample, "FQ1"].values[0],
        fq2=lambda wildcards: samples.loc[samples.Sample==wildcards.sample, "FQ2"].values[0]
    output:
        r1=os.path.join(OUTDIR, "{sample}_R1.trim.fastq.gz"),
        r2=os.path.join(OUTDIR, "{sample}_R2.trim.fastq.gz")
    threads: 4
    shell:
        """
        cutadapt -q 0 -O 1 -m 0 -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
            -o {output.r1} -p {output.r2} {input.fq1} {input.fq2}
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

