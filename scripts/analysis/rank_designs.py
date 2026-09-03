#!/usr/bin/env python3
"""Stage 07 - rank designs.

Reproduces final_ranking.csv from the three cloud co-folder result tables:
  boltz_ipsae_scdockq.csv   (cand, ipsae_min, sc_dockq)
  esm_full_ipsae_scdockq.csv(cand, ipsae_min, sc_dockq)
  protenix_results.csv      (cand, iptm, ..., sc_dockq)   [Protenix has no PAE -> no ipSAE]

pose_PASS = all three sc_DockQ >= 0.23 (dataset threshold).
Among the pose-passing designs, composite z = 4*z(ipSAE_avg) + 1*z(sc_DockQ_avg),
the dataset's 4:1 weighting; ipSAE_avg = mean(Boltz, ESM), sc_DockQ_avg = mean(all 3).

Usage:
  rank_designs.py --scoring DIR [--selection final_selection.csv] [--out final_ranking.csv] [--verify]
"""
import csv, argparse, statistics as st, os, sys

GATE = 0.23

def load(path, keymap):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[r["cand"]] = {k: float(r[v]) for k, v in keymap.items() if r.get(v) not in (None, "")}
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scoring", required=True, help="dir with the 3 co-fold CSVs + final_selection.csv")
    ap.add_argument("--selection", default="final_selection.csv")
    ap.add_argument("--out", default="final_ranking.csv")
    ap.add_argument("--verify", action="store_true", help="compare against an existing --out")
    a = ap.parse_args()
    D = a.scoring

    bz = load(os.path.join(D, "boltz_ipsae_scdockq.csv"),    {"ipsae": "ipsae_min", "scdq": "sc_dockq"})
    ef = load(os.path.join(D, "esm_full_ipsae_scdockq.csv"), {"ipsae": "ipsae_min", "scdq": "sc_dockq"})
    px = load(os.path.join(D, "protenix_results.csv"),       {"scdq": "sc_dockq", "iptm": "iptm"})
    sel = {}
    with open(os.path.join(D, a.selection)) as fh:
        for r in csv.DictReader(fh):
            sel[r["cand"]] = (r.get("backbone", ""), r.get("method", ""))

    cands = [c for c in bz if c in ef and c in px]
    rows = []
    for c in cands:
        bi, ei = bz[c]["ipsae"], ef[c]["ipsae"]
        bs, es, ps = bz[c]["scdq"], ef[c]["scdq"], px[c]["scdq"]
        ipsae_avg = (bi + ei) / 2
        scdq_avg = (bs + es + ps) / 3
        pose = int(bs >= GATE and es >= GATE and ps >= GATE)
        rows.append(dict(cand=c, pose_pass3=pose, ipsae_avg=ipsae_avg, scdq_avg=scdq_avg,
                         bz_ipsae=bi, ef_ipsae=ei, bz_scdq=bs, ef_scdq=es, px_scdq=ps,
                         px_iptm=px[c].get("iptm", ""),
                         backbone=sel.get(c, ("", ""))[0], method=sel.get(c, ("", ""))[1]))

    passed = [r for r in rows if r["pose_pass3"] == 1]
    mi = st.mean(r["ipsae_avg"] for r in passed); si = st.pstdev(r["ipsae_avg"] for r in passed)
    ms = st.mean(r["scdq_avg"] for r in passed);  ss = st.pstdev(r["scdq_avg"] for r in passed)
    for r in passed:
        r["composite"] = 4 * (r["ipsae_avg"] - mi) / si + 1 * (r["scdq_avg"] - ms) / ss
    passed.sort(key=lambda r: -r["composite"])
    for i, r in enumerate(passed, 1):
        r["rank"] = i
    removed = [r for r in rows if r["pose_pass3"] == 0]

    cols = ["rank", "cand", "composite", "pose_pass3", "ipsae_avg", "scdq_avg",
            "bz_ipsae", "ef_ipsae", "bz_scdq", "ef_scdq", "px_scdq", "px_iptm", "backbone", "method"]
    outp = os.path.join(D, a.out)
    if a.verify:
        want = {r["cand"]: r for r in csv.DictReader(open(outp))}
        worst = 0.0
        for r in passed:
            w = want.get(r["cand"])
            if w and w["composite"]:
                worst = max(worst, abs(float(w["composite"]) - r["composite"]))
        print(f"VERIFY: {len(passed)} pose-pass, {len(removed)} removed; "
              f"max |composite delta| vs {a.out} = {worst:.4f}")
        print("OK" if worst < 0.02 else "MISMATCH", "(tolerance 0.02)")
        return

    def fmt(r):
        return {c: (f"{r[c]:.4f}".rstrip("0").rstrip(".") if isinstance(r.get(c), float) else r.get(c, ""))
                for c in cols}
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in passed: w.writerow(fmt(r))
        for r in removed: w.writerow(fmt(r))
    print(f"wrote {outp}: {len(passed)} ranked + {len(removed)} removed")

if __name__ == "__main__":
    main()
