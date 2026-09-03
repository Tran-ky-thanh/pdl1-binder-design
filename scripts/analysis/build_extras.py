import os
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build")
os.makedirs(OUTDIR, exist_ok=True)
import json, base64, glob, io
from Bio.PDB import PDBParser, Superimposer, PDBIO, Structure, Model
OUT=OUTDIR
R=PROJECT
SC=R+"/fold_screen/complex_run/scoring"
CH=R+"/fold_screen/complex_run"
data=json.load(open(OUT+"/report_slim.json"))

# ---- figures for section 02 (metric-code reproduction + consensus validation) ----
# colabfold_consensus.png is produced by make_valfig.py; run it before this script.
for key,path in [("fig_comparison", R+"/comparison.png"),
                 ("fig_consensus", OUT+"/colabfold_consensus.png")]:
    data[key]=base64.b64encode(open(path,"rb").read()).decode()
    print("embedded", key, round(len(data[key])/1024), "KB b64")

# ---- novelty overlay: cand00725 vs nearest reference ----
P=PDBParser(QUIET=True)
def cf(cand):
    for i in range(4):
        g=glob.glob(f"{CH}/chunk{i}/out/{cand}_unrelaxed_rank_001_*.pdb")
        if g: return g[0]
ref_pdb=R+"/files/designs/PD-L1/mythos_preview_single_target_pdl1_rank03/designed/designed.pdb"
cm=P.get_structure("c", cf("cand00725"))[0]   # A=target(115), B=cand binder(72)
rm=P.get_structure("r", ref_pdb)[0]            # A=ref binder(100), B=target(2-110)
ct, rt = cm["A"], rm["B"]                       # PD-L1 in each
cca={r.id[1]:r["CA"] for r in ct if "CA" in r}
fix=[];mov=[]
for r in rt:
    if "CA" in r and r.id[1] in cca:
        fix.append(cca[r.id[1]]); mov.append(r["CA"])
sup=Superimposer(); sup.set_atoms(fix,mov); sup.apply(list(rm.get_atoms()))
print("overlay superpose CA matched=%d rmsd=%.2f"%(len(mov), sup.rms))

news=Structure.Structure("ov"); newm=Model.Model(0); news.add(newm)
tA=ct.copy(); tA.id="A"; newm.add(tA)                 # PD-L1 target
bB=cm["B"].copy(); bB.id="B"; newm.add(bB)            # cand00725 binder
rC=rm["A"].copy(); rC.id="C"; newm.add(rC)            # reference binder (aligned)
bio=PDBIO(); bio.set_structure(news); buf=io.StringIO(); bio.save(buf)
data["novelty_overlay"]="\n".join(l for l in buf.getvalue().splitlines() if l.startswith("ATOM"))+"\n"
data["novelty_overlay_meta"]={"cand":"cand00725","ref":"mythos_preview_single_target_pdl1_rank03",
    "id":36.1,"ref_kd_nM":0.64,"rmsd":round(sup.rms,2)}
print("overlay atoms:", data["novelty_overlay"].count("\n"))

# ---- fold/pose overlay: cand01514 vs its nearest reference BY STRUCTURE (TM=0.815) ----
# Unlike the cand00725 pair above (same epitope, different fold), this is essentially the
# closest structural match across the full top-30 x 78-reference TM-align sweep in
# results/tm_vs_anthropic.csv (global max 0.816): same fold AND same binding pose, still
# <30% sequence identity. rank06 is also an experimentally confirmed binder (Twist Kd 3.0 nM
# / Adaptyv Kd 29.5 nM, both_bind) so it doubles as the affinity example.
ref2_pdb=R+"/files/designs/PD-L1/mythos_preview_single_target_pdl1_rank06/designed/designed.pdb"
cm2=P.get_structure("c2", cf("cand01514"))[0]   # A=target(115), B=cand binder(76)
rm2=P.get_structure("r2", ref2_pdb)[0]           # A=ref binder(72), B=target(109)
ct2, rt2 = cm2["A"], rm2["B"]
cca2={r.id[1]:r["CA"] for r in ct2 if "CA" in r}
fix2=[];mov2=[]
for r in rt2:
    if "CA" in r and r.id[1] in cca2:
        fix2.append(cca2[r.id[1]]); mov2.append(r["CA"])
sup2=Superimposer(); sup2.set_atoms(fix2,mov2); sup2.apply(list(rm2.get_atoms()))
print("fold/pose overlay superpose CA matched=%d rmsd=%.2f"%(len(mov2), sup2.rms))

news2=Structure.Structure("ov2"); newm2=Model.Model(0); news2.add(newm2)
tA2=ct2.copy(); tA2.id="A"; newm2.add(tA2)                # PD-L1 target
bB2=cm2["B"].copy(); bB2.id="B"; newm2.add(bB2)           # cand01514 binder
rC2=rm2["A"].copy(); rC2.id="C"; newm2.add(rC2)           # reference binder (aligned)
bio2=PDBIO(); bio2.set_structure(news2); buf2=io.StringIO(); bio2.save(buf2)
data["foldpose_overlay"]="\n".join(l for l in buf2.getvalue().splitlines() if l.startswith("ATOM"))+"\n"
data["foldpose_overlay_meta"]={"cand":"cand01514","ref":"mythos_preview_single_target_pdl1_rank06",
    "tmscore":0.815,"id":28.9,"ref_kd_nM":3.0,"ref_kd_note":"Twist Kd; Adaptyv Kd 29.5 nM, both_bind",
    "rmsd":round(sup2.rms,2)}
print("fold/pose overlay atoms:", data["foldpose_overlay"].count("\n"))

json.dump(data, open(OUT+"/report_slim.json","w"), separators=(",",":"))
import os
print("report_slim.json MB:", round(os.path.getsize(OUT+"/report_slim.json")/1e6,2))
