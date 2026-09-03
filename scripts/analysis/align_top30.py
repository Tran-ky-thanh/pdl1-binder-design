import io, glob, json
from Bio.PDB import PDBParser, Superimposer, PDBIO
OUT="/mnt/c/Users/thanh/AppData/Local/Temp/claude/C--Users-thanh-Documents/9a8e1c59-098d-4338-a1cc-e43c6d31654e/scratchpad"
CHUNKS="/home/thanh/protein_designs/pdl1_bench/fold_screen/complex_run"
parser=PDBParser(QUIET=True)

def find_pdb(cand):
    for ch in range(4):
        g=glob.glob(f"{CHUNKS}/chunk{ch}/out/{cand}_unrelaxed_rank_001_*.pdb")
        if g: return g[0]
    return None

def target_chain(model):
    # target = chain with the most CA-bearing residues (PD-L1 ~110+ vs binder ~60-80)
    best=None;bn=-1
    for ch in model:
        n=sum(1 for r in ch if 'CA' in r)
        if n>bn: bn=n; best=ch.id
    return best

data=json.load(open(OUT+"/report_slim.json"))
top=[d for d in data["designs"] if d["rank"]<=30]
top.sort(key=lambda d:d["rank"])

# reference = rank-1 target chain, keyed by residue number
refs=parser.get_structure("ref", find_pdb(top[0]["cand"]))
refm=refs[0]; rtc=target_chain(refm)
ref_ca={r.id[1]:r['CA'] for r in refm[rtc] if 'CA' in r}
print("reference:",top[0]["cand"],"target chain",rtc,"CA res",len(ref_ca))

done=0
for d in top:
    p=find_pdb(d["cand"])
    if not p: print("  MISSING",d["cand"]); continue
    s=parser.get_structure(d["cand"], p); m=s[0]; tc=target_chain(m)
    fix=[];mov=[]
    for r in m[tc]:
        if 'CA' in r and r.id[1] in ref_ca:
            fix.append(ref_ca[r.id[1]]); mov.append(r['CA'])
    if len(mov)<20: print("  too few match",d["cand"],len(mov)); continue
    sup=Superimposer(); sup.set_atoms(fix,mov)
    sup.apply(list(m.get_atoms()))   # move WHOLE complex onto the PD-L1 frame
    bio=PDBIO(); bio.set_structure(s); buf=io.StringIO(); bio.save(buf)
    atoms=[l for l in buf.getvalue().splitlines() if l.startswith("ATOM")]
    d["pdb"]="\n".join(atoms)+"\n"
    done+=1
    print(f"  #{d['rank']:>2} {d['cand']} rmsd={sup.rms:.2f} matched={len(mov)} chain={tc}")

# record which chain is the target so the viewer frames it (same length rule client-side is fine)
data["target_chain_ref"]=rtc
json.dump(data, open(OUT+"/report_slim.json","w"), separators=(",",":"))
import os
print("aligned",done,"/",len(top),"  report_slim.json MB:",round(os.path.getsize(OUT+"/report_slim.json")/1e6,2))
