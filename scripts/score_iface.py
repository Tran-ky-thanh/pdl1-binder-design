#!/usr/bin/env python3
"""Orthogonal CPU screen for folded PD-L1 complexes (NO GPU).

Per candidate computes:
  sc_dockq  : DockQ of AF2 complex vs its RFdiffusion designed-backbone pose
              = pose self-consistency (does AF2 dock the binder where the design
              intended?). Native is built from the design backbone with:
                - target chain A renumbered to construct numbering (1-115)
                  (offset 17 for 5C3T/4ZQK/4Z18 author numbering; 0 for pdl1_binder)
                - binder chain B residues RENAMED to the AF2 binder sequence,
                  so DockQ maps chains A:A , B:B cleanly.
              Design binder is poly-Gly -> DockQ / iRMSD / LRMSD are the reliable
              pose signals; fnat is a backbone-only approximation.
  bsa       : buried surface area at interface (Bio.PDB ShrakeRupley), A^2
  n_int_res : interface residues (any cross-chain atom pair < 5 A)
  n_contact : cross-chain residue-residue contacts (< 5 A)

Deps: Biopython + DockQ (already installed in pdl1_bench/.venv). No PyRosetta needed.
Usage: python3 score_iface.py CAND MODEL_PDB
       -> prints one CSV row: cand,sc_dockq,irms,Lrms,fnat,bsa,n_int_res,n_contact
NOTE: the DockQ CLI output parsing (dockq() below) MUST be validated on 2-3
      complexes once tools are reachable -- run with DEBUG=1 to dump raw output.
"""
import os, sys, csv, re, subprocess, tempfile, copy, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

FS = "/home/thanh/protein_designs/pdl1_bench"
SELMAP = FS + "/fold_screen/selection_map.csv"
BBDIRS = [FS + "/rfdiff_out2", FS + "/rfdiff_out"]
P = PDBParser(QUIET=True)

def bb_for(cand):
    with open(SELMAP) as f:
        for row in csv.DictReader(f):
            if row["cand"] == cand:
                return row["backbone"]
    return None

def bb_path(name):
    for d in BBDIRS:
        p = os.path.join(d, name + ".pdb")
        if os.path.exists(p):
            return p
    return None

def load(pdb):
    return P.get_structure("s", pdb)[0]

def binder_resnames(model):
    return {res.id[1]: res.resname for res in model["B"] if res.id[0] == " "}

def build_native(bb_pdb, bbname, bnames, out_pdb):
    m = load(bb_pdb)
    off = 0 if bbname.startswith("pdl1_binder") else 17   # -> construct numbering
    for res in list(m["A"]):
        i = res.id
        res.id = (i[0], i[1] - off, i[2])
    for res in m["B"]:
        if res.id[1] in bnames:
            res.resname = bnames[res.id[1]]
    io = PDBIO(); io.set_structure(m); io.save(out_pdb)

DOCKQ = os.path.join(os.path.dirname(sys.executable), "DockQ")
if not os.path.exists(DOCKQ):
    DOCKQ = "DockQ"

def dockq(model_pdb, native_pdb):
    try:
        r = subprocess.run([DOCKQ, model_pdb, native_pdb, "--mapping", "AB:AB"],
                           capture_output=True, text=True, timeout=180)
    except Exception as e:
        if os.environ.get("DEBUG"):
            sys.stderr.write("DockQ EXC: %r\n" % e)
        return {}
    t = r.stdout + "\n" + r.stderr
    if os.environ.get("DEBUG"):
        sys.stderr.write("=== DockQ raw ===\n" + t + "\n")
    def g(*pats):
        for p in pats:
            m = re.search(p, t)
            if m:
                return float(m.group(1))
        return ""
    return dict(dockq=g(r"DockQ:?\s+([0-9.]+)"),
                irms=g(r"iRMSD:?\s+([0-9.]+)", r"iRMS:?\s+([0-9.]+)"),
                Lrms=g(r"LRMSD:?\s+([0-9.]+)", r"LRMS:?\s+([0-9.]+)"),
                fnat=g(r"[Ff]nat:?\s+([0-9.]+)"))

_sr = ShrakeRupley()
def _sasa_sum(struct):
    _sr.compute(struct, level="A")
    return sum(a.sasa for a in struct.get_atoms())

def bsa(model):
    full = _sasa_sum(copy.deepcopy(model))
    a = copy.deepcopy(model); [a.detach_child(c.id) for c in list(a) if c.id != "A"]
    b = copy.deepcopy(model); [b.detach_child(c.id) for c in list(b) if c.id != "B"]
    return (_sasa_sum(a) + _sasa_sum(b) - full) / 2.0

def contacts(model, cut=5.0):
    ns = NeighborSearch(list(model.get_atoms()))
    pairs, ires = set(), set()
    for a1, a2 in ns.search_all(cut, level="A"):
        r1, r2 = a1.get_parent(), a2.get_parent()
        c1, c2 = r1.get_parent().id, r2.get_parent().id
        if c1 != c2:
            pairs.add((r1.get_full_id(), r2.get_full_id()))
            ires.add(r1.get_full_id()); ires.add(r2.get_full_id())
    return len(ires), len(pairs)

def main():
    cand, model_pdb = sys.argv[1], sys.argv[2]
    sc = {}
    bbname = bb_for(cand)
    bbp = bb_path(bbname) if bbname else None
    if bbp:
        nat = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        try:
            build_native(bbp, bbname, binder_resnames(load(model_pdb)), nat)
            sc = dockq(model_pdb, nat)
        finally:
            try: os.unlink(nat)
            except OSError: pass
    nres, ncon = contacts(load(model_pdb))
    b = bsa(load(model_pdb))
    print("%s,%s,%s,%s,%s,%.1f,%d,%d" % (
        cand, sc.get("dockq", ""), sc.get("irms", ""), sc.get("Lrms", ""),
        sc.get("fnat", ""), b, nres, ncon))

if __name__ == "__main__":
    main()
