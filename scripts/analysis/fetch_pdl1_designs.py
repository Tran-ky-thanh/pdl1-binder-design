#!/usr/bin/env python3
"""Download the FULL Anthropic PD-L1 reference set (`designed` backbone structure for every
design that has one) from the Anthropic/claude-protein-binder-design HuggingFace dataset.

design_summary.parquet lists 90 PD-L1 rows (30 Mythos Preview / multi_target, 30 Mythos
Preview / single_target, 30 Opus 4.8 / multi_target). The current dataset layout keeps each
design's released backbone at:
    data/designs/PD-L1/<full_name>/insilico/designed.cif
Not every row has this file: 12 of the 90 only carry `predicted_*_1to1.cif` (co-fold
predictions against PD-L1), with no released "as-designed" backbone -- those are skipped, so
the practical maximum is 78/90. This script fetches whatever `designed.cif` files are still
missing, converts each to legacy PDB (needed by TM-align / Biopython PDBParser elsewhere in
this repo) and writes it to files/designs/PD-L1/<full_name>/designed/designed.pdb -- the path
scripts/analysis/structural_novelty.py and scripts/analysis/build_extras.py expect.

Usage:
  python scripts/analysis/fetch_pdl1_designs.py
"""
import os
import pandas as pd
from huggingface_hub import hf_hub_download
from Bio.PDB import MMCIFParser, PDBIO

PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
REPO = "Anthropic/claude-protein-binder-design"


def main():
    df = pd.read_parquet(os.path.join(PROJECT, "meta/design_summary.parquet"))
    names = sorted(df[df.target == "PD-L1"]["full_name"].astype(str).tolist())
    print("PD-L1 designs listed in design_summary.parquet:", len(names))

    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    have = fetched = failed = 0

    for n in names:
        local_target = os.path.join(PROJECT, "files/designs/PD-L1", n, "designed", "designed.pdb")
        if os.path.exists(local_target):
            have += 1
            continue
        rel = f"data/designs/PD-L1/{n}/insilico/designed.cif"
        try:
            cif_path = hf_hub_download(REPO, rel, repo_type="dataset")
            struct = parser.get_structure(n, cif_path)
            os.makedirs(os.path.dirname(local_target), exist_ok=True)
            io.set_structure(struct)
            io.save(local_target)
            fetched += 1
            print("got", n)
        except Exception as e:
            failed += 1
            print("no released backbone for", n, "-", type(e).__name__)

    print(f"\nalready had: {have}  fetched: {fetched}  no-backbone (skipped): {failed}  "
          f"total on disk now: {have + fetched} / {len(names)}")


if __name__ == "__main__":
    main()
