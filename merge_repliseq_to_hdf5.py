#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
import h5py
from pyfaidx import Fasta

# -------- Arguments --------
# usage: python merge_repliseq_to_hdf5.py samples.tsv genome.fa outdir hdf5_name
samples_file = sys.argv[1]  # e.g., "samples.tsv"
genome_file = sys.argv[2]   # e.g., "genome.fa"
outdir = sys.argv[3]        # e.g., "results"
hdf5_name = sys.argv[4]     # e.g., "HMEC_repliseq.h5"

os.makedirs(outdir, exist_ok=True)

# -------- Load samples --------
samples = pd.read_csv(samples_file, sep="\t")

# -------- Load genome sizes --------
genome = Fasta(genome_file)
chrom_sizes = {name: len(seq) for name, seq in genome.items()}

# -------- Initialize dense arrays --------
fraction_map = {"Late": 0, "Early": 1, "G1": 2}
chrom_arrays = {chrom: np.zeros((size+1, 3), dtype=np.float32) for chrom, size in chrom_sizes.items()}

# -------- Fill arrays --------
for _, row in samples.iterrows():
    sample, sample_type = row["Sample"], row["Type"]
    cov_file = os.path.join(outdir, f"{sample}.perbase.cov")
    frac_idx = fraction_map[sample_type]

    df = pd.read_csv(cov_file, sep="\t", header=None, names=["chr","pos","depth"])
    df["norm"] = (df["depth"] / df["depth"].sum()) * 4_000_000

    for chrom, subdf in df.groupby("chr"):
        chrom_arrays[chrom][subdf["pos"].values, frac_idx] = subdf["norm"].values

# -------- Write HDF5 --------
hdf5_path = os.path.join(outdir, hdf5_name)
with h5py.File(hdf5_path, "w") as h5:
    for chrom, arr in chrom_arrays.items():
        h5.create_dataset(chrom, data=arr, compression="gzip", chunks=(10000,3))

print(f"[INFO] Wrote {hdf5_path}")
