#!/usr/bin/env python3.11
import argparse, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--specimen_loci", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

df = pd.read_csv(args.specimen_loci, sep="\t")
if "specimen_id" not in df.columns or "locus" not in df.columns:
    df = pd.read_csv(args.specimen_loci, sep="\t", header=None, names=["specimen_id","locus"])

df["v"]=1
wide = (df.pivot_table(index="specimen_id", columns="locus", values="v", aggfunc="max", fill_value=0)
          .reset_index())
wide.to_csv(args.out, sep="\t", index=False)
print(f"[loci_table] wrote {args.out} with {len(wide)} rows")
