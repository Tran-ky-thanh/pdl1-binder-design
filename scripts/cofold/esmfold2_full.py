#!/usr/bin/env python3
"""Stage 06 (ESMFold2-Full) - fold PD-L1:binder complexes via Biohub Forge.

Free, PAE-capable. TWO PASSES with different numpy majors (they cannot coexist):

  # pass 1 - FOLD  (esm SDK needs numpy 2.x)
  pip install -q numpy==2.5.2
  esmfold2_full.py fold  --in final_selection.csv --outdir /tmp/esm_full_struct \
                         --results scoring/esmfold2_full_results.csv [--token ...]

  # pass 2 - SCORE (DockQ needs numpy<2)
  pip install -q numpy==1.26.4
  esmfold2_full.py score --outdir /tmp/esm_full_struct \
                         --out scoring/esm_full_ipsae_scdockq.csv --scoring <scoring_dir>

Token resolution: --token, else $BIOHUB_TOKEN, else /tmp/.bhkey.
'fold' saves <outdir>/<cand>.pdb + <cand>.json (plddt,pae,ptm,iptm) and appends
(cand,esm_full_ptm,esm_full_iptm) to --results. 'score' computes ipSAE (from the json)
and sc_DockQ (from the pdb via score_iface.py) -> cand,ipsae_min,sc_dockq. Both resumable.
"""
import csv, json, os, sys, time, argparse, subprocess, warnings
warnings.filterwarnings("ignore")

DEFAULT_TARGET = ("AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRAR"
                  "LLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNA")

def numpy_major():
    import numpy as np
    return int(np.__version__.split(".")[0])

def get_token(cli):
    if cli: return cli
    if os.environ.get("BIOHUB_TOKEN"): return os.environ["BIOHUB_TOKEN"]
    if os.path.exists("/tmp/.bhkey"): return open("/tmp/.bhkey").read().strip()
    sys.exit("no Biohub token: pass --token, set BIOHUB_TOKEN, or put it in /tmp/.bkey")

def do_fold(a):
    if numpy_major() < 2:
        sys.exit("fold needs numpy>=2 for the esm SDK. Run:  pip install -q numpy==2.5.2")
    import numpy as np
    from esm.sdk.forge import SequenceStructureForgeInferenceClient
    from esm.sdk.api import FoldingConfig, ESMProteinError
    target = open(a.target[1:]).read().split()[0] if a.target.startswith("@") else a.target
    os.makedirs(a.outdir, exist_ok=True)
    rows = [r for r in csv.DictReader(open(a.inp)) if r.get("seq")]
    if not os.path.exists(a.results):
        open(a.results, "w").write("cand,esm_full_ptm,esm_full_iptm\n")
    client = SequenceStructureForgeInferenceClient(model="esmfold2-2026-05", url="https://biohub.ai",
                                                   token=get_token(a.token))
    cfg = FoldingConfig(include_pae=True, include_pair_chains_iptm=True)
    tonp = lambda x: x.detach().cpu().numpy() if hasattr(x, "detach") else np.asarray(x)
    n = 0
    for r in rows:
        c = r["cand"]
        if os.path.exists(os.path.join(a.outdir, c + ".json")):
            continue
        for att in range(3):
            try:
                res = client.fold(target + "|" + r["seq"], config=cfg)
                if isinstance(res, ESMProteinError):
                    print("  ERR %s %r" % (c, res)); break
                pae = tonp(res.pae); pl = tonp(res.plddt); pl = pl[~np.isnan(pl)]
                open(os.path.join(a.outdir, c + ".pdb"), "w").write(res.to_pdb_string())
                json.dump({"plddt": (pl * 100).tolist(), "predicted_aligned_error": pae.tolist(),
                           "ptm": float(res.ptm), "iptm": float(res.interface_ptm)},
                          open(os.path.join(a.outdir, c + ".json"), "w"))
                open(a.results, "a").write("%s,%.4f,%.4f\n" % (c, float(res.ptm), float(res.interface_ptm)))
                n += 1; print("  folded %s iptm=%.3f (%d)" % (c, float(res.interface_ptm), n)); break
            except Exception as e:
                print("  EXC %s att%d %r" % (c, att, e)); time.sleep(5)
        time.sleep(a.space)
    print("fold done, new=%d -> %s" % (n, a.outdir))

def do_score(a):
    if numpy_major() >= 2:
        sys.exit("score needs numpy<2 for DockQ. Run:  pip install -q numpy==1.26.4")
    scripts = os.path.abspath(a.scripts)
    ipsae_py = os.path.join(scripts, "ipsae.py")
    score_py = os.path.join(scripts, "score_iface.py")
    py = sys.executable
    done = {r["cand"] for r in csv.DictReader(open(a.out))} if os.path.exists(a.out) else set()
    if not os.path.exists(a.out): open(a.out, "w").write("cand,ipsae_min,sc_dockq\n")
    cands = sorted(f[:-5] for f in os.listdir(a.outdir) if f.endswith(".json"))
    for c in cands:
        if c in done: continue
        pdb = os.path.join(a.outdir, c + ".pdb"); js = os.path.join(a.outdir, c + ".json")
        subprocess.run([py, ipsae_py, js, pdb, "10", "10"], capture_output=True, text=True, timeout=180)
        txt = os.path.join(a.outdir, c + "_10_10.txt"); ip = ""
        if os.path.exists(txt):
            vals = [float(p[5]) for p in (l.split() for l in open(txt)) if len(p) > 5 and p[4] == "asym"]
            ip = "%.4f" % min(vals) if vals else ""
        try:
            out = subprocess.check_output([py, score_py, c, pdb], stderr=subprocess.DEVNULL, timeout=300).decode()
            row = [l for l in out.splitlines() if l.startswith(c + ",")]
            sq = row[-1].split(",")[1] if row else ""
        except Exception as e:
            print("  sc_dockq failed %s: %r" % (c, e)); sq = ""
        open(a.out, "a").write("%s,%s,%s\n" % (c, ip, sq))
        print("  scored %s ipSAE=%s sc_dockq=%s" % (c, ip, sq))
    print("score done -> %s" % a.out)

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fold"); f.set_defaults(fn=do_fold)
    f.add_argument("--in", dest="inp", required=True)
    f.add_argument("--outdir", default="/tmp/esm_full_struct")
    f.add_argument("--results", required=True)
    f.add_argument("--token", default="")
    f.add_argument("--target", default=DEFAULT_TARGET)
    f.add_argument("--space", type=float, default=0.4)
    s = sub.add_parser("score"); s.set_defaults(fn=do_score)
    s.add_argument("--outdir", default="/tmp/esm_full_struct")
    s.add_argument("--out", required=True)
    s.add_argument("--scoring", required=True)
    s.add_argument("--scripts", default=os.path.dirname(os.path.abspath(__file__)) + "/..")
    a = ap.parse_args()
    a.fn(a)

if __name__ == "__main__":
    main()
