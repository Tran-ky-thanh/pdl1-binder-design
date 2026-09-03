import os
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build")
os.makedirs(OUTDIR, exist_ok=True)
import json, csv
OUT=OUTDIR
SC=PROJECT+"/fold_screen/complex_run/scoring"
data=json.load(open(OUT+"/report_data.json"))
# strip pdb from ranks >30 to keep file lean
for d in data["designs"]:
    if d["rank"]>30: d["pdb"]=None
# all 48 rows (incl 5 removed) for pose-pass scatter
all48=[]
for r in csv.DictReader(open(SC+"/final_ranking.csv")):
    all48.append(dict(cand=r["cand"], pose=int(r["pose_pass3"]),
        composite=(float(r["composite"]) if r["composite"] else None),
        bz_ipsae=float(r["bz_ipsae"]), ef_ipsae=float(r["ef_ipsae"]),
        bz_scdq=float(r["bz_scdq"]), ef_scdq=float(r["ef_scdq"]),
        px_scdq=float(r["px_scdq"]), px_iptm=float(r["px_iptm"])))
data["all48"]=all48
json.dump(data, open(OUT+"/report_slim.json","w"), separators=(",",":"))
import os
print("report_slim.json MB:", round(os.path.getsize(OUT+"/report_slim.json")/1e6,2))
print("top30 with pdb:", sum(1 for d in data["designs"] if d["pdb"]))
print("all48 rows:", len(all48), "removed:", sum(1 for r in all48 if r["pose"]==0))
# quick: nearest-neighbour confirmed-binder note for top30
nb_conf=sum(1 for d in data["designs"][:30] if d["nov_isbinder"])
print("top30 whose nearest dataset design is a confirmed binder:", nb_conf)
print("max nov_id top30:", max(d["nov_id"] for d in data["designs"][:30]))
print("max nov_id all43:", max(d["nov_id"] for d in data["designs"]))
