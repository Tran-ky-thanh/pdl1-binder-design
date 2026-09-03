#!/usr/bin/env python3
"""Stage 05 (aggregate) — interface metrics for every co-folded complex.

Production driver for the 1,050 co-folds. For each complex (chain A = PD-L1 target,
chain B = binder) it computes, under numpy 1.26.4 (PyRosetta + DockQ need numpy<2):

  ipSAE_min                      via ipsae.py, from the PAE json (pure-PAE)
  sc_DockQ, BSA, contacts        via score_iface.py (DockQ vs the designed backbone)
  ΔΔG, dSASA, CMS, n_int_res     via PyRosetta InterfaceAnalyzer

and writes one aggregated table in the ranked_1050.csv / iface_1050.csv schema, which the
two-stage screen (filter_stages.py: Stage-1 ipSAE+sc_DockQ, Stage-2 ΔΔG/CMS physics) consumes.

The PyRosetta ΔΔG/CMS core (constrained full relax + InterfaceAnalyzer + CMS) is the SAME
one used by the QC benchmark `bench_ddg.py`; this script imports it from there rather than
duplicating it. Lesson (why the *constrained full* relax): recycle-1 AF2 structures carry
small backbone clashes; a side-chain-only relax leaves them in place and yields false
"clash blow-ups" in ΔΔG. A coordinate-constrained FastRelax with the backbone free removes
the artefact — ipSAE/sc_DockQ barely move, ΔΔG/CMS only become meaningful after it.

Resumable (one .row file per design). Parallelise with --slice K N (process cands[K::N]) and
run N copies, then merge with --aggregate-only.

Usage:
  export PDL1_PROJECT=/path/to/working/data          # holds fold_screen/complex_run/chunk*/out/*.pdb
  python scripts/score_complexes.py --out build/iface_1050.csv            # all, sequential
  for k in 0 1 2 3 4 5; do python scripts/score_complexes.py --slice $k 6 & done; wait
  python scripts/score_complexes.py --aggregate-only --out build/iface_1050.csv
"""
import os, sys, glob, csv, subprocess, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)          # so we can reuse bench_ddg.py's PyRosetta ddG/CMS core
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build"); os.makedirs(OUTDIR, exist_ok=True)
IPSAE = os.path.join(HERE, "ipsae.py")
SCORE = os.path.join(HERE, "score_iface.py")
SUFF = "_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000"
COLS = ["cand", "ipsae_min", "sc_dockq", "bsa", "ddg", "dSASA", "cms", "nres", "irms", "Lrms", "fnat", "ncontact"]

def _f(x):
    try: float(x); return True
    except (TypeError, ValueError): return False

def ipsae_min(pdb, pae, py):
    txt = pdb[:-4] + "_10_10.txt"
    if not os.path.exists(txt):
        subprocess.run([py, IPSAE, pae, pdb, "10", "10"], capture_output=True, text=True, timeout=300)
    if not os.path.exists(txt): return ""
    vals = [float(p[5]) for p in (l.split() for l in open(txt)) if len(p) > 5 and p[4] == "asym" and _f(p[5])]
    return "%.4f" % min(vals) if vals else ""

def sc_iface(cand, pdb, py):
    try:
        r = subprocess.run([py, SCORE, cand, pdb], capture_output=True, text=True, timeout=300)
        p = r.stdout.strip().split(",")
        return p[1], p[2], p[3], p[4], p[5], p[6], p[7]   # scq,irms,Lrms,fnat,bsa,n_int,n_contact
    except Exception:
        return ("",) * 7

def make_scorer():
    import bench_ddg as bd           # module-level pyrosetta.init(...) + relax()/cms()/InterfaceAnalyzerMover
    def score(pdb, rlx_out):
        pose = bd.pose_from_pdb(pdb)
        bd.relax(pose)               # coordinate-constrained full relax (the ddG lesson)
        pose.dump_pdb(rlx_out)
        sf = bd.pyrosetta.get_fa_scorefxn()
        ia = bd.InterfaceAnalyzerMover(1, False, sf, False, False, True)   # jump 1 = A|B interface
        ia.apply(pose)
        cm = bd.cms(pose)
        return ("%.2f" % ia.get_separated_interface_energy(),
                "%.1f" % ia.get_interface_delta_sasa(),
                ("%.1f" % cm if cm != "" else ""),
                str(ia.get_num_interface_residues()))
    return score

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default=os.path.join(PROJECT, "fold_screen/complex_run"),
                    help="dir with chunk*/out/<cand>...pdb complexes")
    ap.add_argument("--out", default=os.path.join(OUTDIR, "iface_1050.csv"))
    ap.add_argument("--rowdir", default=os.path.join(OUTDIR, "rows_iface"))
    ap.add_argument("--slice", nargs=2, type=int, metavar=("K", "N"), help="process cands[K::N]")
    ap.add_argument("--aggregate-only", action="store_true", help="just merge existing .row files")
    a = ap.parse_args()
    os.makedirs(a.rowdir, exist_ok=True)
    rlxdir = os.path.join(OUTDIR, "relaxed"); os.makedirs(rlxdir, exist_ok=True)

    pdbs = sorted(glob.glob(os.path.join(a.indir, "chunk*/out/*%s.pdb" % SUFF)))
    cands = sorted(os.path.basename(p).split("_unrelaxed")[0] for p in pdbs)

    if not a.aggregate_only:
        py = sys.executable
        mine = cands[a.slice[0]::a.slice[1]] if a.slice else cands
        score = make_scorer()
        for cand in mine:
            rowf = os.path.join(a.rowdir, cand + ".row")
            if os.path.exists(rowf): continue
            g = glob.glob(os.path.join(a.indir, "chunk*/out/%s%s.pdb" % (cand, SUFF)))
            if not g: continue
            pdb = g[0]; pae = pdb.replace(SUFF + ".pdb", "_predicted_aligned_error_v1.json")
            try:
                ip = ipsae_min(pdb, pae, py)
                rlx = os.path.join(rlxdir, cand + ".pdb")
                ddg, dsasa, cm, nres = score(pdb, rlx)
                scq, irms, lrms, fnat, bsa, nir, ncon = sc_iface(cand, rlx, py)
                tmp = rowf + ".tmp"
                with open(tmp, "w") as f:
                    f.write(",".join([cand, ip, scq, bsa, ddg, dsasa, cm, nres, irms, lrms, fnat, ncon]) + "\n")
                os.rename(tmp, rowf)
                print("scored", cand, "ddg=%s cms=%s ipSAE=%s scDockQ=%s" % (ddg, cm, ip, scq))
            except Exception as e:
                sys.stderr.write("ERR %s: %r\n" % (cand, e))

    rows = []
    for cand in cands:
        rf = os.path.join(a.rowdir, cand + ".row")
        if os.path.exists(rf):
            rows.append(next(csv.reader(open(rf))))
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(COLS); w.writerows(rows)
    print("wrote %s (%d designs)" % (a.out, len(rows)))

if __name__ == "__main__":
    main()
