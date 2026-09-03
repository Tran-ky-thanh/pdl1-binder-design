import json, base64, glob, io
from Bio.PDB import PDBParser, Superimposer, PDBIO, Structure, Model
OUT="/mnt/c/Users/thanh/AppData/Local/Temp/claude/C--Users-thanh-Documents/9a8e1c59-098d-4338-a1cc-e43c6d31654e/scratchpad"
R="/home/thanh/protein_designs/pdl1_bench"
SC=R+"/fold_screen/complex_run/scoring"
CH=R+"/fold_screen/complex_run"
data=json.load(open(OUT+"/report_slim.json"))

# ---- figures for section 02 ----
for key,path in [("fig_comparison", R+"/comparison.png"), ("fig_relax", SC+"/relax_compare.png")]:
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

json.dump(data, open(OUT+"/report_slim.json","w"), separators=(",",":"))
import os
print("report_slim.json MB:", round(os.path.getsize(OUT+"/report_slim.json")/1e6,2))
