#!/usr/bin/env python
"""Merge kallisto abundance.tsv files into TPM and estimated-count matrices."""
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--samples", nargs="+", required=True)
p.add_argument("--inputs", nargs="+", required=True)
p.add_argument("--tpm", required=True)
p.add_argument("--counts", required=True)
a = p.parse_args()
assert len(a.samples) == len(a.inputs)

tpm, cnt = [], []
for s, f in zip(a.samples, a.inputs):
    df = pd.read_csv(f, sep="\t", index_col="target_id")
    tpm.append(df["tpm"].rename(s))
    cnt.append(df["est_counts"].rename(s))
for series, out in ((tpm, a.tpm), (cnt, a.counts)):
    m = pd.concat(series, axis=1)
    m.index.name = "transcript_id"
    m.to_csv(out, sep="\t")
print("done")
