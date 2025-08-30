#!/usr/bin/env python3.11
import argparse
from pathlib import Path
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Pivot specimen×locus table; keep flexible column names.")
    ap.add_argument("--specimen_loci", required=True, help="TSV with specimen/locus/seq_id/length")
    ap.add_argument("--out", required=True, help="Output TSV pivot (presence/absence or max length)")
    ap.add_argument("--taxon_map", default="data/staging/taxon_map.tsv",
                    help="Specimen->pretty species map (will be created/updated if missing)")
    args = ap.parse_args()

    df = pd.read_csv(args.specimen_loci, sep="\t", dtype=str)

    # Normalize column names
    cols = {c.lower(): c for c in df.columns}
    if "specimen_id" in cols:
        spec_col = cols["specimen_id"]
    elif "specimen" in cols:
        spec_col = cols["specimen"]
        df = df.rename(columns={spec_col: "specimen_id"})
        spec_col = "specimen_id"
    else:
        raise SystemExit("specimen_loci TSV must have a 'specimen' or 'specimen_id' column.")

    if "locus" not in df.columns:
        # try lowercase fallbacks
        if "Locus" in df.columns: df = df.rename(columns={"Locus": "locus"})
        else: raise SystemExit("specimen_loci TSV must have a 'locus' column.")

    # presence=1 per specimen×locus (or use length if you prefer)
    df["v"] = 1
    piv = (df
           .pivot_table(index="specimen_id", columns="locus", values="v",
                        fill_value=0, aggfunc="max")
           .reset_index())

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    piv.to_csv(args.out, sep="\t", index=False)

    # Ensure/refresh a minimal taxon_map (only add missing specimen_ids as 'Unknown')
    tmap_path = Path(args.taxon_map)
    tmap_path.parent.mkdir(parents=True, exist_ok=True)
    if tmap_path.exists():
        tmap = pd.read_csv(tmap_path, sep="\t", header=None, names=["specimen_id","species"], dtype=str)
    else:
        tmap = pd.DataFrame(columns=["specimen_id","species"], dtype=str)

    known = set(tmap["specimen_id"].astype(str))
    missing = [s for s in piv["specimen_id"].astype(str) if s not in known]
    if missing:
        add = pd.DataFrame({"specimen_id": missing, "species": ["Unknown"]*len(missing)})
        tmap = pd.concat([tmap, add], ignore_index=True)
        tmap.drop_duplicates(subset=["specimen_id"], keep="first", inplace=True)
        tmap.to_csv(tmap_path, sep="\t", header=False, index=False)
        print(f"[enrich] taxon_map: added {len(missing)} new specimen ids to {tmap_path}")
    else:
        print("[enrich] taxon_map: no new specimen ids to add")

if __name__ == "__main__":
    main()
