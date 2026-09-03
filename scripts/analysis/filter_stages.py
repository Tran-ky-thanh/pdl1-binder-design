#!/usr/bin/env python3
"""Two-stage screen of the 1,050 co-folded designs (Stage-1 -> Stage-2).

Stage-1 (1,050 -> 242): ipSAE + sc_DockQ consensus gate; its survivors are recorded in
                        stage1_pass.txt.
Stage-2 (242 -> 157):   PyRosetta *physics* gate applied on top of Stage-1, using the
                        interface energetics measured by PyRosetta InterfaceAnalyzer for
                        all 1,050 co-folds (columns `ddg`, `cms` in ranked_1050.csv):
                        keep designs with  ΔΔG <= -40 REU  AND  CMS >= 360 Å².

This removes 85 designs that pass the co-fold consensus but have weak/repulsive interfaces
(a cheap physical-plausibility layer before the expensive cloud co-folding stage). The
survivors then go to curation (liability + clustering -> 50).

Runs from a clean clone against the committed tables:
  python scripts/analysis/filter_stages.py --ranked results/ranked_1050.csv \
         --stage1 results/stage1_pass.txt --verify
Thresholds are overridable: --ddg -40 --cms 360
"""
import csv, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ranked", default="results/ranked_1050.csv",
                    help="per-design interface metrics incl. ddg/cms for all 1,050")
    ap.add_argument("--stage1", default="results/stage1_pass.txt",
                    help="Stage-1 survivors (ipSAE + sc_DockQ)")
    ap.add_argument("--out", default="results/stage2_pass.txt")
    ap.add_argument("--ddg", type=float, default=-40.0, help="max ΔΔG (REU)")
    ap.add_argument("--cms", type=float, default=360.0, help="min contact molecular surface (Å²)")
    ap.add_argument("--verify", action="store_true", help="compare against an existing --out")
    a = ap.parse_args()

    R = {r["cand"]: r for r in csv.DictReader(open(a.ranked))}
    s1 = [x.strip() for x in open(a.stage1) if x.strip()]
    def fnum(c, k):
        try: return float(R[c][k])
        except (KeyError, ValueError, TypeError): return None

    passed = [c for c in s1 if c in R and (fnum(c, "ddg") or 0) <= a.ddg and (fnum(c, "cms") or 0) >= a.cms]
    print("Stage-1 in: %d  ->  Stage-2 (ΔΔG <= %.0f REU & CMS >= %.0f Å²): %d  (PyRosetta physics removed %d)"
          % (len(s1), a.ddg, a.cms, len(passed), len(s1) - len(passed)))

    if a.verify:
        want = set(x.strip() for x in open(a.out) if x.strip())
        got = set(passed)
        print("VERIFY vs %s: exact match = %s  (extra=%d, missing=%d)"
              % (a.out, got == want, len(got - want), len(want - got)))
        return
    with open(a.out, "w") as f:
        for c in passed:
            f.write(c + "\n")
    print("wrote %s (%d designs)" % (a.out, len(passed)))

if __name__ == "__main__":
    main()
