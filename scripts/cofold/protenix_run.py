#!/usr/bin/env python3
"""Stage 06 (Protenix) - fold PD-L1:binder complexes with Protenix v2 via JapanFold.

Free public API (no key), but RATE-LIMITED, so we submit in small batches (5-8 jobs),
space submits a few seconds apart, wait for the batch to finish, then cool down 2-3 min
before the next batch. Resumable: candidates already in --out are skipped.

Protenix exposes no PAE -> no ipSAE. We record iptm/ptm/plddt/confidence from results.json
and compute sc_DockQ from the predicted CIF via score_iface.py (DockQ vs designed backbone).

Output CSV columns: cand,iptm,ptm,plddt,confidence,sc_dockq

Usage:
  protenix_run.py --in final_selection.csv --out scoring/protenix_results.csv \
                  --scoring <scoring_dir> [--batch 6] [--space 3.5] [--cooldown 150]
"""
import requests, csv, time, io, zipfile, json, os, sys, argparse, subprocess, warnings
warnings.filterwarnings("ignore")

BASE = "https://api.japanfold.aiand.com"
HDRS = {  # JapanFold sits behind Cloudflare - browser-like headers are required
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36",
    "Accept": "application/json", "Content-Type": "application/json",
    "Referer": "https://japanfold.aiand.com/",
}
# PD-L1 IgV target construct (chain A) used for every co-fold in this campaign.
DEFAULT_TARGET = ("AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRAR"
                  "LLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNA")

def submit(name, target, binder, retries=4):
    inp = ("sequences:\n  - protein: {id: A, sequence: %s}\n"
           "  - protein: {id: B, sequence: %s}\n" % (target, binder))
    body = {"model": "protenix-v2", "name": name, "input": inp}
    for k in range(retries):
        try:
            r = requests.post(BASE + "/v1/predictions", headers=HDRS, json=body, timeout=60)
            if r.status_code in (200, 201, 202):
                return r.json().get("id")
            # 429 / Cloudflare rate-limit -> exponential backoff and retry
            print("    submit HTTP %s (retry %d): %s" % (r.status_code, k + 1, r.text[:120]))
        except Exception as e:
            print("    submit error (retry %d): %r" % (k + 1, e))
        time.sleep(8 * (k + 1))
    return None

def poll(jid):
    try:
        s = requests.get(BASE + "/v1/jobs/" + jid, headers=HDRS, timeout=40).json()
        return s.get("status", "?"), bool(s.get("done") or s.get("finished_at"))
    except Exception:
        return "err", False

def fetch(jid, out_cif):
    z = requests.get(BASE + "/v1/jobs/" + jid + "/archive", headers=HDRS, timeout=120).content
    with zipfile.ZipFile(io.BytesIO(z)) as zf:
        names = zf.namelist()
        res = json.loads(zf.read([n for n in names if n.endswith("results.json")][0]))[0]
        cif = [n for n in names if n.endswith(".cif")][0]
        with open(out_cif, "wb") as fh:
            fh.write(zf.read(cif))
    return res

def cif_to_pdb(cif, pdb):
    from Bio.PDB import MMCIFParser, PDBIO
    s = MMCIFParser(QUIET=True).get_structure("s", cif)
    io_ = PDBIO(); io_.set_structure(s); io_.save(pdb)

def sc_dockq(cand, pdb, scripts_dir):
    """Run score_iface.py (DockQ vs the designed backbone). Returns sc_dockq or ''."""
    try:
        out = subprocess.check_output([sys.executable, os.path.join(scripts_dir, "score_iface.py"), cand, pdb],
                                      stderr=subprocess.DEVNULL, timeout=180).decode().strip().splitlines()
        # last CSV line: cand,sc_dockq,irms,Lrms,fnat,bsa,n_int_res,n_contact
        row = [l for l in out if l.startswith(cand + ",")]
        return row[-1].split(",")[1] if row else ""
    except Exception as e:
        print("    sc_dockq failed for %s: %r" % (cand, e))
        return ""

def chunks(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="CSV with columns cand,seq")
    ap.add_argument("--out", required=True, help="results CSV (appended, resumable)")
    ap.add_argument("--scoring", required=True, help="dir holding score_iface.py's deps context")
    ap.add_argument("--scripts", default=os.path.dirname(os.path.abspath(__file__)) + "/..",
                    help="scripts/ dir containing score_iface.py")
    ap.add_argument("--outdir", default="/tmp/protenix_cif", help="where to save predicted CIF/PDB")
    ap.add_argument("--target", default=DEFAULT_TARGET, help="target sequence, or @file")
    ap.add_argument("--batch", type=int, default=6, help="jobs submitted per batch (5-8)")
    ap.add_argument("--space", type=float, default=3.5, help="seconds between submits (avoid rapid-submit block)")
    ap.add_argument("--cooldown", type=float, default=150, help="seconds to wait between batches (2-3 min)")
    ap.add_argument("--batch-timeout", type=float, default=600, help="max seconds to wait for a batch")
    a = ap.parse_args()

    target = open(a.target[1:]).read().split()[0] if a.target.startswith("@") else a.target
    scripts = os.path.abspath(a.scripts)
    os.makedirs(a.outdir, exist_ok=True)

    done = set()
    if os.path.exists(a.out):
        with open(a.out) as fh:
            done = {r["cand"] for r in csv.DictReader(fh)}
    todo = [r for r in csv.DictReader(open(a.inp)) if r["cand"] not in done and r.get("seq")]
    print("candidates: %d to run, %d already done" % (len(todo), len(done)))

    new = not os.path.exists(a.out)
    fout = open(a.out, "a", newline="")
    w = csv.writer(fout)
    if new:
        w.writerow(["cand", "iptm", "ptm", "plddt", "confidence", "sc_dockq"]); fout.flush()

    for bi, batch in enumerate(chunks(todo, a.batch), 1):
        print("\n=== batch %d (%d jobs) ===" % (bi, len(batch)))
        jobs = {}
        for r in batch:
            jid = submit(r["cand"], target, r["seq"])
            print("  submit %s -> %s" % (r["cand"], jid))
            if jid:
                jobs[r["cand"]] = jid
            time.sleep(a.space)  # pace submits within the batch
        # wait for this batch to finish
        pending = dict(jobs); t0 = time.time()
        while pending and time.time() - t0 < a.batch_timeout:
            time.sleep(15)
            for cand, jid in list(pending.items()):
                st, fin = poll(jid)
                if fin:
                    cif = os.path.join(a.outdir, cand + ".cif")
                    pdb = os.path.join(a.outdir, cand + ".pdb")
                    try:
                        res = fetch(jid, cif)
                        cif_to_pdb(cif, pdb)
                        scdq = sc_dockq(cand, pdb, scripts)
                        w.writerow([cand, res.get("iptm", ""), res.get("ptm", ""),
                                    res.get("plddt", res.get("complex_plddt", "")),
                                    res.get("confidence_score", res.get("confidence", "")), scdq])
                        fout.flush()
                        print("  done %s  iptm=%s sc_dockq=%s (t+%ds)" %
                              (cand, res.get("iptm"), scdq, int(time.time() - t0)))
                    except Exception as e:
                        print("  fetch/score failed %s: %r" % (cand, e))
                    del pending[cand]
        if pending:
            print("  batch timeout, still pending: %s" % list(pending))
        if bi * a.batch < len(todo):
            print("  cooldown %.0fs before next batch..." % a.cooldown)
            time.sleep(a.cooldown)

    fout.close()
    print("\nwrote %s" % a.out)

if __name__ == "__main__":
    main()
