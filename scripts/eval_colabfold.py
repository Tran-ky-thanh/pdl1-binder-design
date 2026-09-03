#!/usr/bin/env python3
"""
Evaluate localColabFold predictions (--num-recycle 5, no templates) for the 6
selected PD-L1 designs, and compare against the dataset's good/bad picture.

For each design's ColabFold rank_001 model we compute:
  * ColabFold's own confidence: ipTM, pTM, mean pLDDT (whole + per chain)
  * ipSAE (same ipsae.py path, cutoffs 10/10) on the CF model + CF PAE json
  * sc_DockQ of the CF complex vs the released design model (same binder-backbone
    / target-heavy recipe) -> "did CF fold the intended pose?"
Then line them up with the dataset consensus ipsae_min and the wet-lab verdict.
"""
import os
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build")
os.makedirs(OUTDIR, exist_ok=True)
import os, glob, json, sys
import numpy as np
import pandas as pd
sys.path.insert(0, PROJECT)
import recompute_pdl1 as R

WORK = PROJECT
OUT = os.path.join(WORK, "colabfold", "out2")


def cf_files(name):
    pdb = sorted(glob.glob(f"{OUT}/{name}/{name}_unrelaxed_rank_001_*.pdb"))
    js = sorted(glob.glob(f"{OUT}/{name}/{name}_scores_rank_001_*.json"))
    return (pdb[0] if pdb else None), (js[0] if js else None)


def main():
    sel = pd.read_csv(os.path.join(WORK, "colabfold", "selection.csv"))
    pred = pd.read_parquet(os.path.join(WORK, "predictions.parquet"))
    des_map = (pred[pred.target == "PD-L1"]
               .drop_duplicates("design_name")
               .set_index("design_name").designed_path.to_dict())

    recs = []
    for _, s in sel.iterrows():
        name, bl = s.full_name, int(s.binder_length)
        pdb, js = cf_files(name)
        if not pdb or not js:
            print(f"MISSING CF output for {name}")
            continue
        d = json.load(open(js))
        plddt = np.array(d["plddt"], dtype=float)

        szp = R.chain_sizes(pdb)
        bpred = R.pick_binder_chain(szp, bl)
        tpred = [c for c in szp if c != bpred][0]
        # per-chain pLDDT (chain order in pdb = A target, B binder)
        nA = szp.get("A", 0)
        plddt_t = plddt[:nA].mean() if nA else np.nan
        plddt_b = plddt[nA:].mean() if nA else np.nan

        ip = R.run_ipsae(js, pdb)          # AF2 mode (json + pdb)
        b2t = ip[(bpred, tpred)]
        t2b = ip[(tpred, bpred)]
        ipmin = min(b2t, t2b)

        scdockq = np.nan
        dp = des_map.get(name, "")
        if isinstance(dp, str) and dp.endswith(".cif"):
            R.extract_members([dp])
            des_cif = os.path.join(R.FILES, dp)
            des = des_cif[:-4] + ".pdb"
            R.cif_to_pdb(des_cif, des)
            szd = R.chain_sizes(des)
            bdes = R.pick_binder_chain(szd, bl)
            mf = pdb[:-4] + "_binderBB.pdb"
            nf = des[:-4] + "_binderBB.pdb"
            R.write_filtered(pdb, bpred, mf)
            R.write_filtered(des, bdes, nf)
            scdockq = R.run_dockq(mf, nf)["dockq"]

        recs.append(dict(
            full_name=name, tier=s.tier, binder_len=bl,
            ds_consensus=s.consensus, ds_af3of3=s.af3of3, ds_boltz2=s.boltz2,
            wetlab_binder=s.binder_final, kd_nM=s.kd_nM,
            cf_iptm=round(float(d.get("iptm", np.nan)), 3),
            cf_ptm=round(float(d.get("ptm", np.nan)), 3),
            cf_plddt=round(plddt.mean(), 1),
            cf_plddt_binder=round(plddt_b, 1),
            cf_ipsae_min=round(ipmin, 3),
            cf_ipsae_b2t=round(b2t, 3),
            cf_scdockq_vs_design=round(scdockq, 3),
        ))
        print(f"{s.tier:14s} {name:42s} "
              f"CF ipTM={d.get('iptm'):.3f} ipSAEmin={ipmin:.3f} "
              f"scDockQ={scdockq:.3f}  (ds cons {s.consensus})", flush=True)

    res = pd.DataFrame(recs)
    res.to_csv(os.path.join(WORK, "colabfold", "cf_eval.csv"), index=False)
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)
    print("\n===== ColabFold (num-recycle 5, no templates) vs dataset =====")
    cols = ["tier", "wetlab_binder", "kd_nM", "ds_consensus",
            "cf_iptm", "cf_ipsae_min", "cf_plddt_binder", "cf_scdockq_vs_design"]
    print(res.sort_values("ds_consensus", ascending=False)[cols].to_string(index=False))
    print(f"\nwrote colabfold/cf_eval.csv")


if __name__ == "__main__":
    main()
