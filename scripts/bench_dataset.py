#!/usr/bin/env python3
"""Reference BSA / interface-geometry distribution from the Anthropic dataset's
predicted complexes (default predictor: afm3 = AF-multimer v3, same family as our
AF2 folds -> directly comparable). Chain-agnostic: uses the first two chains.
ddG/CMS is NOT here (needs PyRosetta -> use rosetta_iface.py on these same .cif).

Usage: python3 bench_dataset.py [PREDICTOR]     # default afm3
       -> prints CSV: design,pred,bsa,n_int_res,n_contact,cA_nres,cB_nres
"""
import sys, glob, os, copy, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

D = "/home/thanh/protein_designs/pdl1_bench/files/designs/PD-L1"
PRED = sys.argv[1] if len(sys.argv) > 1 else "afm3"
_sr = ShrakeRupley()

def sasa_sum(s):
    _sr.compute(s, level="A")
    return sum(a.sasa for a in s.get_atoms())

def load(f):
    return MMCIFParser(QUIET=True).get_structure("x", f)[0]

def bsa_generic(model):
    chains = [c.id for c in model]
    if len(chains) < 2:
        return None
    c1, c2 = chains[0], chains[1]
    full = sasa_sum(copy.deepcopy(model))
    a = copy.deepcopy(model); [a.detach_child(c.id) for c in list(a) if c.id != c1]
    b = copy.deepcopy(model); [b.detach_child(c.id) for c in list(b) if c.id != c2]
    return (sasa_sum(a) + sasa_sum(b) - full) / 2.0, c1, c2

def contacts(model, c1, c2, cut=5.0):
    ns = NeighborSearch(list(model.get_atoms()))
    pr, ir = set(), set()
    for a1, a2 in ns.search_all(cut, level="A"):
        x = a1.get_parent().get_parent().id
        y = a2.get_parent().get_parent().id
        if {x, y} == {c1, c2}:
            r1, r2 = a1.get_parent(), a2.get_parent()
            pr.add((r1.get_full_id(), r2.get_full_id()))
            ir.add(r1.get_full_id()); ir.add(r2.get_full_id())
    return len(ir), len(pr)

print("design,pred,bsa,n_int_res,n_contact,cA_nres,cB_nres")
for d in sorted(glob.glob(D + "/*/")):
    name = os.path.basename(d.rstrip("/"))
    cands = glob.glob(os.path.join(d, PRED, "seed_*", "model.cif"))
    if not cands:
        continue
    f = sorted(cands)[0]
    try:
        m = load(f)
        r = bsa_generic(m)
        if not r:
            continue
        b, c1, c2 = r
        nres, ncon = contacts(m, c1, c2)
        ns = {c.id: len([x for x in c if x.id[0] == " "]) for c in m}
        print("%s,%s,%.1f,%d,%d,%d,%d" % (name, PRED, b, nres, ncon, ns.get(c1, 0), ns.get(c2, 0)))
    except Exception as e:
        sys.stderr.write("ERR %s: %r\n" % (name, e))
