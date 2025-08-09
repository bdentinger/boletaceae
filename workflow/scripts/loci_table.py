#!/usr/bin/env python3
"""
Create a specimen × locus presence/absence table from specimen_loci.tsv.
Input:  data/staging/specimens/specimen_loci.tsv
Output: data/staging/loci_table.tsv
"""
import argparse
import pandas as pd
from pathlib import Path

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('--specimen_loci', default='data/staging/specimens/specimen_loci.tsv')
  ap.add_argument('--out', default='data/staging/loci_table.tsv')
  args = ap.parse_args()

  df = pd.read_csv(args.specimen_loci, sep='\t')
  pt = (df.assign(v=1)
          .pivot_table(index='specimen_id', columns='locus', values='v', fill_value=0, aggfunc='max')
          .reset_index())
  Path(args.out).parent.mkdir(parents=True, exist_ok=True)
  pt.to_csv(args.out, sep='\t', index=False)

if __name__ == "__main__":
  main()
