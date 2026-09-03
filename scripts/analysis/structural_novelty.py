#!/usr/bin/env python3
"""Structural novelty of our top-30 binders vs the FULL Anthropic PD-L1 reference set.

Sequence-level novelty (scripts/analysis/novelty_and_extract.py) already compares every
candidate's sequence against all 90 Anthropic PD-L1 rows in design_summary.parquet (the
parquet carries `sequence` for every row, no PDB download needed). This script complements
that with a STRUCTURAL comparison, which does need the actual 3-D coordinates: for each of
our top-30 binders, TM-align its isolated binder chain against every Anthropic PD-L1 binder
chain for which a `designed.pdb` structure is present locally, and keep the best (highest
TM-score) match.

Requires `files/designs/PD-L1/*/designed/designed.pdb` on disk (download via
scripts/analysis/fetch_pdl1_designs.py) and a compiled TMalign binary (TM-align, Zhang Lab,
zhanggroup.org/TM-align — build with `g++ -O2 -o TMalign TMalign.cpp`; not installed by
requirements.txt because it is a standalone C++ binary, not a Python package).

Usage:
  python scripts/analysis/structural_novelty.py \
      --tmalign /home/thanh/bin/TMalign --top 30 \
      --out results/tm_vs_anthropic.csv
"""
import os
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build")
os.makedirs(OUTDIR, exist_ok=True)

import argparse, csv, glob, re, statistics as st, subprocess, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select

AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
      'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
      'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
PDL1_SIG = "VVEYGSNMTIEC"  # PD-L1 IgV-domain signature; whichever chain lacks it is the binder
P = PDBParser(QUIET=True)


def seq(chain):
    return "".join(AA.get(r.resname, "") for r in chain if 'CA' in r)


class OneChain(Select):
    def __init__(self, cid): self.cid = cid
    def accept_chain(self, c): return c.id == self.cid


def save_chain(model, cid, out):
    io = PDBIO(); io.set_structure(model); io.save(out, OneChain(cid))


def find_our_pdb(cand, complex_run):
    for ch in range(4):
        g = glob.glob(f"{complex_run}/chunk{ch}/out/{cand}_unrelaxed_rank_001_*.pdb")
        if g:
            return g[0]
    return None


def tmscore(tmalign_bin, q, t):
    try:
        out = subprocess.check_output([tmalign_bin, q, t], text=True, timeout=60,
                                       stderr=subprocess.DEVNULL)
    except Exception:
        return None
    vals = re.findall(r"TM-score=\s*([0-9.]+)", out)
    return max(float(v) for v in vals) if vals else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmalign", default=os.environ.get("TMALIGN_BIN", "TMalign"),
                     help="path to the compiled TMalign binary")
    ap.add_argument("--ranking", default="results/final_ranking.csv")
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--anthropic-glob",
                     default=os.path.join(PROJECT, "files/designs/PD-L1/*/designed/designed.pdb"))
    ap.add_argument("--out", default="results/tm_vs_anthropic.csv")
    a = ap.parse_args()

    complex_run = os.path.join(PROJECT, "fold_screen/complex_run")
    workdir = os.path.join(OUTDIR, "tmalign")
    ours_dir = os.path.join(workdir, "ours")
    anth_dir = os.path.join(workdir, "anthropic")
    os.makedirs(ours_dir, exist_ok=True)
    os.makedirs(anth_dir, exist_ok=True)

    top = []
    for r in csv.DictReader(open(a.ranking)):
        if r.get("rank") and int(r["rank"]) <= a.top:
            top.append((int(r["rank"]), r["cand"]))
    top.sort()

    for rank, cand in top:
        pdb = find_our_pdb(cand, complex_run)
        if not pdb:
            continue
        m = P.get_structure("c", pdb)[0]
        binder = [c.id for c in m if PDL1_SIG not in seq(c)]
        if binder:
            save_chain(m, binder[0], os.path.join(ours_dir, cand + ".pdb"))
    print("our binder chains extracted:", len(glob.glob(ours_dir + "/*.pdb")))

    anth_paths = sorted(glob.glob(a.anthropic_glob))
    for d in anth_paths:
        name = os.path.basename(os.path.dirname(os.path.dirname(d)))
        m = P.get_structure("a", d)[0]
        binder = [c.id for c in m if PDL1_SIG not in seq(c) and len(seq(c)) > 20]
        if binder:
            save_chain(m, binder[0], os.path.join(anth_dir, name + ".pdb"))
    anth_files = sorted(glob.glob(anth_dir + "/*.pdb"))
    print("Anthropic PD-L1 reference binder chains:", len(anth_files),
          "(full dataset has 90 PD-L1 designs; this uses however many `designed.pdb` "
          "files are present on disk)")

    rows = []
    for rank, cand in top:
        q = os.path.join(ours_dir, cand + ".pdb")
        if not os.path.exists(q):
            continue
        best = (0.0, None)
        for t in anth_files:
            s = tmscore(a.tmalign, q, t)
            if s and s > best[0]:
                best = (s, os.path.basename(t)[:-4])
        rows.append((rank, cand, best[0], best[1]))
        print(f"#{rank:>2} {cand}  best TM={best[0]:.3f} vs {best[1]}")

    tms = [r[2] for r in rows if r[2]]
    if tms:
        print("\n=== structural novelty: top-%d binder fold vs %d Anthropic PD-L1 binders ==="
              % (a.top, len(anth_files)))
        print("TM-score to nearest Anthropic binder: median=%.3f min=%.3f max=%.3f"
              % (st.median(tms), min(tms), max(tms)))
        print(">=0.5 (same fold): %d/%d   0.4-0.5: %d   <0.4 (distinct fold): %d"
              % (sum(t >= 0.5 for t in tms), len(tms),
                 sum(0.4 <= t < 0.5 for t in tms), sum(t < 0.4 for t in tms)))

    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "cand", "best_tmscore", "nearest_anthropic"])
        w.writerows(rows)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
