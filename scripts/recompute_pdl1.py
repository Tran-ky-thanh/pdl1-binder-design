#!/usr/bin/env python3
"""
Recompute ipSAE_min and sc_DockQ for the PD-L1 subset of
Anthropic/claude-protein-binder-design and compare to the released values.

Pipeline reverse-engineered from the dataset's column_dictionary:
  * ipSAE : Dunbrack ipsae.py, the "ipSAE" column (d0res variant),
            pae_cutoff=10, dist_cutoff=10.
            binder_to_target = align on binder, score over target.
            ipsae_min = min(binder_to_target, target_to_binder).
  * sc_DockQ : DockQ v2 of predicted complex (model) vs design model (native),
            with binder = backbone atoms only (N,CA,C,O), target = heavy atoms,
            best chain assignment. -> matches sc_dockq / dockq_fnat/irms/lrms.

Structures + PAE live inside a 74 GB zip in the repo; we pull only the members
we need via HTTP range requests (remotezip). Nothing else is downloaded whole.
"""
import os, re, sys, subprocess, argparse, warnings
import numpy as np
import pandas as pd
import gemmi
warnings.simplefilter("ignore")
from Bio.PDB import PDBParser, PDBIO, Select
from remotezip import RemoteZip
from huggingface_hub import hf_hub_download

REPO = "Anthropic/claude-protein-binder-design"
ZIP_URL = (f"https://huggingface.co/datasets/{REPO}/resolve/main/"
           "structure_and_pae/protein_binder_design_structure_and_pae_release.zip")
ZIP_PREFIX = "protein_binder_design_structure_and_pae_release/"
BB = {"N", "CA", "C", "O"}
WORK = "/home/thanh/protein_designs/pdl1_bench"
FILES = os.path.join(WORK, "files")
os.makedirs(FILES, exist_ok=True)
PDB = PDBParser(QUIET=True)


def extract_members(members):
    """Download the given zip-internal relative paths via range requests (cached)."""
    todo = sorted({m for m in members
                   if m and not os.path.exists(os.path.join(FILES, m))})
    if not todo:
        return
    print(f"  extracting {len(todo)} files from remote zip ...", flush=True)
    with RemoteZip(ZIP_URL) as z:
        for m in todo:
            data = z.read(ZIP_PREFIX + m)
            out = os.path.join(FILES, m)
            os.makedirs(os.path.dirname(out), exist_ok=True)
            with open(out, "wb") as fh:
                fh.write(data)


def cif_to_pdb(cif, pdb):
    """Robust CIF -> PDB via gemmi (handles every model's mmCIF dialect)."""
    st = gemmi.read_structure(cif)
    st.setup_entities()
    with open(pdb, "w") as fh:
        fh.write(st.make_pdb_string())


def chain_sizes(pdb):
    m = PDB.get_structure("x", pdb)[0]
    return {ch.id: sum(1 for _ in ch.get_residues()) for ch in m}


def pick_binder_chain(sizes, binder_len):
    """chain whose residue count == binder_length; else the smaller chain."""
    hits = [c for c, n in sizes.items() if n == binder_len]
    if len(hits) == 1:
        return hits[0]
    return min(sizes, key=sizes.get)


def write_filtered(pdb_in, binder_chain, out):
    """binder chain -> backbone atoms only; other chain(s) -> heavy atoms."""
    m = PDB.get_structure("x", pdb_in)[0]

    class Sel(Select):
        def accept_atom(self, a):
            if a.element == "H":
                return 0
            ch = a.get_parent().get_parent().id
            if ch == binder_chain:
                return 1 if a.get_name() in BB else 0
            return 1

    io = PDBIO()
    io.set_structure(m)
    io.save(out, Sel())


def run_ipsae(pae_path, pdb_path):
    """Return {(chn1,chn2): ipSAE_d0res} for the two asym rows."""
    subprocess.run([sys.executable, "ipsae.py", pae_path, pdb_path, "10", "10"],
                   cwd=WORK, capture_output=True, text=True)
    txt = os.path.splitext(pdb_path)[0] + "_10_10.txt"   # ipsae.py output naming
    out = {}
    with open(txt) as fh:
        for line in fh:
            p = line.split()
            if len(p) > 6 and p[4] == "asym":
                out[(p[0], p[1])] = float(p[5])  # col 5 == "ipSAE" (d0res)
    return out


def run_dockq(model_pdb, native_pdb):
    r = subprocess.run(["DockQ", model_pdb, native_pdb],
                       cwd=WORK, capture_output=True, text=True)
    o = r.stdout

    def g(key):
        m = re.search(rf"{key}:\s*([0-9.]+)", o)
        return float(m.group(1)) if m else np.nan

    tot = re.search(r"Total DockQ.*?:\s*([0-9.]+)", o)
    return dict(dockq=float(tot.group(1)) if tot else np.nan,
                fnat=g("fnat"), irms=g("iRMSD"), lrms=g("LRMSD"))


def build_sample(per_model, seed):
    df = pd.read_parquet(os.path.join(WORK, "predictions.parquet"))
    p = df[(df.target == "PD-L1")
           & (df.designed_path.astype(str).str.endswith(".cif"))].copy()
    # binder_length via design_summary joined on uuid
    ds_path = hf_hub_download(REPO, "data/tables/design_summary.parquet",
                              repo_type="dataset")
    ds = pd.read_parquet(ds_path)[["uuid", "binder_length"]]
    p = p.merge(ds, on="uuid", how="left")
    # one representative row per (model, design): prefer requested seed, else first
    p["seed_str"] = p.seed.astype(str)
    p["_pref"] = (p.seed_str == str(seed)).astype(int)
    p = p.sort_values(["cofolding_model", "design_name", "_pref", "seed_str"],
                      ascending=[True, True, False, True])
    reps = p.drop_duplicates(["cofolding_model", "design_name"], keep="first")
    parts = [g.sort_values("design_name").head(per_model)
             for _, g in reps.groupby("cofolding_model")]
    return pd.concat(parts).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-model", type=int, default=2)
    ap.add_argument("--seed", default="1")
    ap.add_argument("--out", default=os.path.join(WORK, "comparison.csv"))
    args = ap.parse_args()

    sample = build_sample(args.per_model, args.seed)
    print(f"sample rows: {len(sample)}  "
          f"(models={sample.cofolding_model.nunique()}, "
          f"designs={sample.design_name.nunique()})", flush=True)

    need = set()
    for _, r in sample.iterrows():
        need |= {r.path_to_structure, r.path_to_pae, r.designed_path}
    extract_members(need)

    recs = []
    for i, r in sample.iterrows():
        struct_cif = os.path.join(FILES, r.path_to_structure)
        pae = os.path.join(FILES, r.path_to_pae)
        des_cif = os.path.join(FILES, r.designed_path)
        try:
            struct = struct_cif[:-4] + ".pdb"
            des = des_cif[:-4] + ".pdb"
            cif_to_pdb(struct_cif, struct)
            cif_to_pdb(des_cif, des)

            szp = chain_sizes(struct)
            szd = chain_sizes(des)
            bl = int(r.binder_length)
            bpred = pick_binder_chain(szp, bl)
            bdes = pick_binder_chain(szd, bl)
            tpred = [c for c in szp if c != bpred][0]

            ip = run_ipsae(pae, struct)
            b2t = ip[(bpred, tpred)]              # align on binder
            t2b = ip[(tpred, bpred)]              # align on target
            ipsae_min = min(b2t, t2b)

            mf = struct[:-4] + "_binderBB.pdb"
            nf = des[:-4] + "_binderBB.pdb"
            write_filtered(struct, bpred, mf)
            write_filtered(des, bdes, nf)
            dq = run_dockq(mf, nf)

            recs.append(dict(
                design_name=r.design_name, cofolding_model=r.cofolding_model,
                seed=r.seed, binder_len=bl,
                my_ipsae_min=ipsae_min, ref_ipsae_min=float(r.ipsae_min),
                my_b2t=b2t, ref_b2t=float(r.ipsae_binder_to_target),
                my_t2b=t2b, ref_t2b=float(r.ipsae_target_to_binder),
                my_scdockq=dq["dockq"], ref_scdockq=float(r.sc_dockq),
                my_fnat=dq["fnat"], ref_fnat=float(r.dockq_fnat),
                my_irms=dq["irms"], ref_irms=float(r.dockq_irms),
                my_lrms=dq["lrms"], ref_lrms=float(r.dockq_lrms),
            ))
            print(f"[{i+1}/{len(sample)}] {r.cofolding_model:8s} {r.design_name[:34]:34s} "
                  f"ipSAEmin {ipsae_min:.3f}/{r.ipsae_min:.3f}  "
                  f"scDockQ {dq['dockq']:.3f}/{r.sc_dockq:.3f}", flush=True)
        except Exception as e:
            print(f"[{i+1}] FAIL {r.design_name} {r.cofolding_model}: {e}", flush=True)

    res = pd.DataFrame(recs)
    res.to_csv(args.out, index=False)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(11, 5))
        for ax, (a, b, name) in zip(
                axs, [("my_ipsae_min", "ref_ipsae_min", "ipSAE_min"),
                      ("my_scdockq", "ref_scdockq", "sc_DockQ")]):
            for m, g in res.groupby("cofolding_model"):
                ax.scatter(g[b], g[a], s=18, label=m, alpha=0.8)
            lim = [min(res[a].min(), res[b].min()) - 0.02, 1.02]
            ax.plot(lim, lim, "k--", lw=0.8)
            ax.set_xlim(lim); ax.set_ylim(lim)
            ax.set_xlabel(f"released {name}")
            ax.set_ylabel(f"recomputed {name}")
            r = res[a].corr(res[b]); md = (res[a] - res[b]).abs().max()
            ax.set_title(f"{name}   r={r:.5f}   max|Δ|={md:.4f}")
        axs[1].legend(fontsize=7, ncol=2, loc="lower right", title="cofolding model")
        fig.suptitle(f"PD-L1 subset: recomputed vs released "
                     f"(n={len(res)}, {res.cofolding_model.nunique()} cofolding models)")
        fig.tight_layout()
        fig.savefig(os.path.join(WORK, "comparison.png"), dpi=140)
        print("wrote comparison.png")
    except Exception as e:
        print("plot skipped:", e)

    print("\n===== AGREEMENT (mine vs released) =====")
    for a, b, name in [("my_ipsae_min", "ref_ipsae_min", "ipsae_min"),
                       ("my_scdockq", "ref_scdockq", "sc_dockq"),
                       ("my_fnat", "ref_fnat", "dockq_fnat"),
                       ("my_irms", "ref_irms", "dockq_irms"),
                       ("my_lrms", "ref_lrms", "dockq_lrms")]:
        d = (res[a] - res[b]).abs()
        r2 = res[a].corr(res[b])
        print(f"{name:12s} n={len(res):3d}  max|Δ|={d.max():.4f}  "
              f"mean|Δ|={d.mean():.4f}  pearson_r={r2:.5f}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
