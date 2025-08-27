#!/usr/bin/env python3.11
import argparse, re, sys
from pathlib import Path

def parse_header(h: str):
    """
    Expect headers like:  BD619-Boletus-cervinococcineus
    Returns:
      specimen_key = BD619
      species_pretty = Boletus cervinococcineus
      species_label  = Boletus-cervinococcineus
    """
    h = h.strip()
    if h.startswith(">"):
        h = h[1:]
    # split on first hyphen into key and species-ish
    if "-" in h:
        key, rest = h.split("-", 1)
    else:
        # fallback: whole header is species-ish, synthesize a key
        key, rest = h, h
    # species label: keep hyphens as-is
    species_label = rest.strip()
    # pretty species (spaces)
    species_pretty = species_label.replace("-", " ").strip()
    specimen_key = key.strip()
    return specimen_key, species_pretty, species_label

ap = argparse.ArgumentParser(description="Ingest custom FASTA to internal format and update taxon_map.")
ap.add_argument("--in", dest="inp", required=True, help="input custom FASTA")
ap.add_argument("--locus", required=True, help="RPB1|RPB2|TEF1|ITS|LSU ...")
ap.add_argument("--out", required=True, help="output FASTA in internal header format or preserved")
ap.add_argument("--taxon-map-append", required=True, help="append mapping to this TSV (specimen\\tspecies)")
ap.add_argument("--preserve-labels", action="store_true",
                help="write FASTA headers exactly as the original label (to match constraint tips).")
args = ap.parse_args()

INP = Path(args.inp)
OUT = Path(args.out)
MAP = Path(args.taxon_map_append)
OUT.parent.mkdir(parents=True, exist_ok=True)
MAP.parent.mkdir(parents=True, exist_ok=True)

written = 0
with INP.open("r", encoding="utf-8", errors="ignore") as fh_in, \
     OUT.open("w", encoding="utf-8") as fh_out, \
     MAP.open("a", encoding="utf-8") as fh_map:

    hdr = None
    seq = []
    for line in fh_in:
        if line.startswith(">"):
            # flush previous
            if hdr is not None and seq:
                key, sp_pretty, sp_label = parse_header(hdr)
                if args.preserve_labels:
                    # keep exactly the species-like label as header
                    new_id = sp_label
                else:
                    # internal format: specimen|LOCUS|Custom:<key>
                    new_id = f"{key}|{args.locus}|Custom:{key}"
                fh_out.write(f">{new_id}\n{''.join(seq)}")
                # taxon map: specimen -> pretty species
                fh_map.write(f"{key}\t{sp_pretty}\n")
                written += 1
            hdr = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    # flush last
    if hdr is not None and seq:
        key, sp_pretty, sp_label = parse_header(hdr)
        if args.preserve_labels:
            new_id = sp_label
        else:
            new_id = f"{key}|{args.locus}|Custom:{key}"
        fh_out.write(f">{new_id}\n{''.join(seq)}")
        fh_map.write(f"{key}\t{sp_pretty}\n")
        written += 1

print(f"[ingest_custom_fasta] wrote {written} records to {OUT}", file=sys.stderr)
