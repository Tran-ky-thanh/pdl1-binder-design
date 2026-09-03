#!/usr/bin/env python3
"""Stage 06 (Boltz-2) - fold PD-L1:binder complexes with Boltz-2 via api.boltz.bio.

Needs a LIVE key (sk_bc_ws_live_...); test keys (sk_bc_ws_test_) return a MOCK fixture.
Submit -> poll -> download the prediction archive -> compute ipSAE (from sample_0_pae.npz)
and sc_DockQ (from the CIF via score_iface.py). Resumable via the jobs file.

Runs entirely under numpy 1.26.4 (Boltz is an HTTP API; ipSAE + DockQ need numpy<2).

Output CSV columns: cand,ipsae_min,sc_dockq   (plus a metrics file with iptm/ptm/plddt).

Usage:
  boltz2_run.py --in final_selection.csv --out scoring/boltz_ipsae_scdockq.csv \
                --scoring <scoring_dir> [--key sk_bc_ws_live_...] [--space 0.7]
Key resolution: --key, else $BOLTZ_API_KEY, else /tmp/.bkey.
"""
import requests, csv, time, io, tarfile, os, sys, argparse, subprocess, warnings
warnings.filterwarnings("ignore")

BASE = "https://api.boltz.bio/compute/v1/predictions/structure-and-binding"
DEFAULT_TARGET = ("AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRAR"
                  "LLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNA")

def get_key(cli):
    if cli: return cli
    if os.environ.get("BOLTZ_API_KEY"): return os.environ["BOLTZ_API_KEY"]
    if os.path.exists("/tmp/.bkey"): return open("/tmp/.bkey").read().strip()
    sys.exit("no Boltz key: pass --key, set BOLTZ_API_KEY, or put it in /tmp/.bkey")

def body(target, binder):
    return {"input": {"entities": [
        {"chain_ids": ["A"], "type": "protein", "value": target},
        {"chain_ids": ["B"], "type": "protein", "value": binder, "msa": {"type": "empty"}}],
        "binding": {"type": "protein_protein_binding", "binder_chain_ids": ["B"]},
        "num_samples": 1}, "model": "boltz-2.1"}

def api(H, url, method="GET", js=None, tries=4):
    for i in range(tries):
        try:
            r = requests.request(method, url, headers=H, json=js, timeout=90)
            if r.status_code in (200, 201):
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(3 * (i + 1)); continue
            print("  API %s %s: %s" % (method, r.status_code, r.text[:150])); return None
        except Exception as e:
            print("  API err (%d): %r" % (i, e)); time.sleep(3)
    return None

def cif_to_pdb(cif, pdb):
    from Bio.PDB import MMCIFParser, PDBIO
    s = MMCIFParser(QUIET=True).get_structure("s", cif)
    io_ = PDBIO(); io_.set_structure(s); io_.save(pdb)

def ipsae_min(py, ipsae_py, npz, cif):
    subprocess.run([py, ipsae_py, npz, cif, "10", "10"], capture_output=True, text=True, timeout=180)
    txt = cif[:-4] + "_10_10.txt"
    if os.path.exists(txt):
        vals = [float(p[5]) for p in (l.split() for l in open(txt)) if len(p) > 5 and p[4] == "asym"]
        if vals: return "%.4f" % min(vals)
    return ""

def sc_dockq(py, score_py, cand, pdb):
    try:
        out = subprocess.check_output([py, score_py, cand, pdb], stderr=subprocess.DEVNULL, timeout=300).decode()
        row = [l for l in out.splitlines() if l.startswith(cand + ",")]
        return row[-1].split(",")[1] if row else ""
    except Exception as e:
        print("    sc_dockq failed %s: %r" % (cand, e)); return ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="CSV with cand,seq")
    ap.add_argument("--out", required=True, help="cand,ipsae_min,sc_dockq (resumable)")
    ap.add_argument("--scoring", required=True, help="dir with ipsae.py/score_iface.py context")
    ap.add_argument("--scripts", default=os.path.dirname(os.path.abspath(__file__)) + "/..")
    ap.add_argument("--jobs", default="", help="jobs csv (cand,jobid); default <out dir>/boltz_jobs.csv")
    ap.add_argument("--metrics", default="", help="metrics csv; default <out dir>/boltz2_results.csv")
    ap.add_argument("--outdir", default="/tmp/boltz_struct")
    ap.add_argument("--target", default=DEFAULT_TARGET)
    ap.add_argument("--key", default="")
    ap.add_argument("--space", type=float, default=0.7, help="seconds between submits")
    ap.add_argument("--poll", type=float, default=20, help="seconds between poll cycles")
    ap.add_argument("--timeout", type=float, default=1800, help="max seconds to wait for all jobs")
    a = ap.parse_args()

    H = {"x-api-key": get_key(a.key), "Content-Type": "application/json"}
    scripts = os.path.abspath(a.scripts)
    ipsae_py = os.path.join(scripts, "ipsae.py")
    score_py = os.path.join(scripts, "score_iface.py")
    target = open(a.target[1:]).read().split()[0] if a.target.startswith("@") else a.target
    os.makedirs(a.outdir, exist_ok=True)
    outdir_csv = os.path.dirname(os.path.abspath(a.out))
    jobs_f = a.jobs or os.path.join(outdir_csv, "boltz_jobs.csv")
    metr_f = a.metrics or os.path.join(outdir_csv, "boltz2_results.csv")
    py = sys.executable

    rows = [r for r in csv.DictReader(open(a.inp)) if r.get("seq")]
    jobmap = {r["cand"]: r["jobid"] for r in csv.DictReader(open(jobs_f))} if os.path.exists(jobs_f) else {}
    scored = {r["cand"] for r in csv.DictReader(open(a.out))} if os.path.exists(a.out) else set()

    # 1) submit missing
    if not os.path.exists(jobs_f): open(jobs_f, "w").write("cand,jobid\n")
    for r in rows:
        c = r["cand"]
        if c in jobmap: continue
        res = api(H, BASE, "POST", body(target, r["seq"]))
        if res and res.get("id"):
            jobmap[c] = res["id"]
            open(jobs_f, "a").write("%s,%s\n" % (c, res["id"]))
            print("  submit %s -> %s" % (c, res["id"]))
        else:
            print("  SUBMIT FAIL %s" % c)
        time.sleep(a.space)

    # 2) poll + 3) score as each finishes
    if not os.path.exists(a.out): open(a.out, "w").write("cand,ipsae_min,sc_dockq\n")
    if not os.path.exists(metr_f):
        open(metr_f, "w").write("cand,status,iptm,ptm,complex_plddt,structure_confidence\n")
    pending = {c: j for c, j in jobmap.items() if c not in scored and c in [r["cand"] for r in rows]}
    t0 = time.time()
    while pending and time.time() - t0 < a.timeout:
        for c, jid in list(pending.items()):
            st = api(H, BASE + "/" + jid)
            if not st: continue
            status = st.get("status")
            if status not in ("succeeded", "failed", "cancelled", "error"):
                continue
            if status != "succeeded":
                print("  %s %s" % (c, status)); del pending[c]; continue
            m = (st.get("output") or {}).get("best_sample", {}).get("metrics", {})
            open(metr_f, "a").write("%s,%s,%s,%s,%s,%s\n" % (
                c, status, m.get("iptm", ""), m.get("ptm", ""),
                m.get("complex_plddt", ""), m.get("structure_confidence", "")))
            try:
                dd = os.path.join(a.outdir, c); os.makedirs(dd, exist_ok=True)
                arch = st["output"]["archive"]["url"]
                tgz = os.path.join(dd, "a.tar.gz")
                open(tgz, "wb").write(requests.get(arch, timeout=180).content)
                with tarfile.open(tgz) as t: t.extractall(dd)
                cif = os.path.join(dd, "prediction", "sample_0_predicted_structure.cif")
                npz = os.path.join(dd, "prediction", "sample_0_pae.npz")
                ip = ipsae_min(py, ipsae_py, npz, cif)
                pdb = os.path.join(dd, "model.pdb"); cif_to_pdb(cif, pdb)
                sq = sc_dockq(py, score_py, c, pdb)
                open(a.out, "a").write("%s,%s,%s\n" % (c, ip, sq))
                print("  done %s  ipSAE=%s sc_dockq=%s  iptm=%s" % (c, ip, sq, m.get("iptm", "")))
            except Exception as e:
                print("  score failed %s: %r" % (c, e))
            del pending[c]
        if pending:
            print("  ...waiting on %d jobs (t+%ds)" % (len(pending), int(time.time() - t0)))
            time.sleep(a.poll)
    if pending: print("  timeout, still pending: %s" % list(pending))
    print("wrote %s" % a.out)

if __name__ == "__main__":
    main()
