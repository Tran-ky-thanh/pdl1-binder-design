import os
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build")
os.makedirs(OUTDIR, exist_ok=True)
import csv, os, glob, json
import pandas as pd
SC=PROJECT+"/fold_screen/complex_run/scoring"
OUT=OUTDIR
CHUNKS=PROJECT+"/fold_screen/complex_run"

# our 43 ranked
rank={}
for r in csv.DictReader(open(SC+"/final_ranking.csv")):
    if r.get("rank"): rank[r["cand"]]=r
sel={r["cand"]:r for r in csv.DictReader(open(SC+"/final_selection.csv"))}

# anthropic PD-L1 binders
df=pd.read_parquet(PROJECT+"/meta/design_summary.parquet")
print("targets:",df['target'].value_counts().to_dict())
pdl1=df[df['target'].astype(str).str.contains('PD-L1',case=False,na=False)].copy()
print("PD-L1 rows:",len(pdl1))
conf=pdl1[pdl1['binder_final']==True]
print("PD-L1 confirmed binders:",len(conf),"kd range nM:",
      round(conf['kd_nM_final'].min(),1) if len(conf) else None,
      round(conf['kd_nM_final'].max(),1) if len(conf) else None)
anth=[(row.full_name,str(row.sequence),bool(row.binder_final),row.kd_nM_final,row.rank)
      for row in pdl1.itertuples() if isinstance(row.sequence,str) and len(row.sequence)>10]
print("anth PD-L1 seqs:",len(anth))

from Bio.Align import PairwiseAligner, substitution_matrices
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load("BLOSUM62")
al.open_gap_score=-11; al.extend_gap_score=-1; al.mode="global"
def ident(q,t):
    try: aln=al.align(q,t)[0]
    except Exception: return 0.0
    s=str(aln).split("\n"); # count matches from aligned rows
    a=aln[0]; b=aln[1]; m=sum(1 for x,y in zip(a,b) if x==y and x!='-')
    return 100.0*m/len(q)

def find_pdb(cand):
    for ch in range(4):
        g=glob.glob(f"{CHUNKS}/chunk{ch}/out/{cand}_unrelaxed_rank_001_*.pdb")
        if g: return g[0]
    return None

def atomonly(path):
    out=[]
    for ln in open(path):
        if ln.startswith("ATOM"):
            out.append(ln.rstrip("\n"))
    return "\n".join(out)+"\n"

designs=[]
for cand,r in sorted(rank.items(), key=lambda kv:int(kv[1]["rank"])):
    seq=sel[cand]["seq"]
    # novelty
    best=(0.0,None,None,None)
    for fn,ts,bf,kd,rk in anth:
        i=ident(seq,ts)
        if i>best[0]: best=(i,fn,bf,kd)
    pdb=find_pdb(cand)
    designs.append(dict(
        rank=int(r["rank"]), cand=cand, composite=float(r["composite"]),
        ipsae_avg=float(r["ipsae_avg"]), scdq_avg=float(r["scdq_avg"]),
        bz_ipsae=float(r["bz_ipsae"]), ef_ipsae=float(r["ef_ipsae"]),
        bz_scdq=float(r["bz_scdq"]), ef_scdq=float(r["ef_scdq"]), px_scdq=float(r["px_scdq"]),
        backbone=r["backbone"], method=r["method"], seq=seq, length=len(seq),
        nov_id=round(best[0],1), nov_match=best[1], nov_isbinder=best[2],
        nov_kd=(round(best[3],1) if best[3]==best[3] and best[3] is not None else None),
        pdb=(atomonly(pdb) if pdb else None)
    ))
    print(f"#{r['rank']:>2} {cand} id={best[0]:5.1f}% match={best[1]} pdb={'Y' if pdb else 'N'}")

data=dict(
    funnel=dict(cofold=1050, stage1=242, stage2=157, selected=50, threemodel=48, posepass=43),
    anth_pdl1=len(pdl1), anth_conf=len(conf),
    anth_kd_min=(round(float(conf['kd_nM_final'].min()),1) if len(conf) else None),
    anth_kd_max=(round(float(conf['kd_nM_final'].max()),1) if len(conf) else None),
    designs=designs)
json.dump(data, open(OUT+"/report_data.json","w"))
print("\nwrote report_data.json  size MB:", round(os.path.getsize(OUT+"/report_data.json")/1e6,2))
print("with pdb:", sum(1 for d in designs if d['pdb']), "/", len(designs))
