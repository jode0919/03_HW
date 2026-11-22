#!/usr/bin/env python
"""Merge per-sample featureCounts tables into one gene x sample matrix."""
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--samples", nargs="+", required=True)
p.add_argument("--inputs", nargs="+", required=True)
p.add_argument("--out", required=True)
a = p.parse_args()
assert len(a.samples) == len(a.inputs)

cols = []
for s, f in zip(a.samples, a.inputs):
    df = pd.read_csv(f, sep="\t", comment="#", index_col="Geneid")
    cols.append(df.iloc[:, -1].rename(s))   # last column = counts
mat = pd.concat(cols, axis=1)
mat.index.name = "gene_id"
mat.to_csv(a.out, sep="\t")
print(f"Wrote {mat.shape[0]} genes x {mat.shape[1]} samples")
