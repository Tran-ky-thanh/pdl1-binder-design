#!/usr/bin/env python3
"""ddG / CMS / SC via PyRosetta InterfaceAnalyzer (Baker minibinder metrics).
Runs on dataset afm3 (.cif, converted to PDB) + a random sample of our AF2 (.pdb).
Interface A_B (both sets: chain A=target 115, chain B=binder). FastRelax(1) then
analyze (RELAX=0 to skip). Writes rows incrementally (flush) so partial results
are readable via the file mount while it runs.
Usage: python3 bench_ddg.py OUT_CSV [N_OURS=20]
"""
import os, sys, glob, tempfile, random, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, PDBIO
import pyrosetta
pyrosetta.init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

FS = "/home/thanh/protein_designs/pdl1_bench"
DATA = FS + "/files/designs/PD-L1"
RUN = FS + "/fold_screen/complex_run"

def cif_to_pdb(cif):
    s = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    p = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    io = PDBIO(); io.set_structure(s); io.save(p)
    return p

def relax(pose):
    from pyrosetta.rosetta.protocols.relax import FastRelax
    sf = pyrosetta.get_fa_scorefxn()
    fr = FastRelax(sf, 1)
    try: fr.constrain_relax_to_start_coords(True)  # keep backbone near AF2 pose
    except Exception: pass
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(True); mm.set_chi(True)  # full relax (bb+chi), constrained to start
    fr.set_movemap(mm); fr.apply(pose)

def cms(pose):
    try:
        S = pyrosetta.rosetta.core.select.residue_selector.ChainSelector
        F = pyrosetta.rosetta.protocols.simple_filters.ContactMolecularSurfaceFilter
        f = F(); f.selector1(S("A")); f.selector2(S("B"))
        return round(f.compute(pose), 1)
    except Exception:
        return ""

def score(path):
    tmp = None
    if path.endswith(".cif"):
        path = tmp = cif_to_pdb(path)
    try:
        pose = pose_from_pdb(path)
        if os.environ.get("RELAX", "1") == "1":
            relax(pose)
        sf = pyrosetta.get_fa_scorefxn()
        # jump 1 = interface between the 2 chains (A=target, B=binder);
        # args: (interface_jump, tracer, sf, compute_packstat, pack_input, pack_separated)
        ia = InterfaceAnalyzerMover(1, False, sf, False, False, True)
        ia.apply(pose)
        ddg = ia.get_separated_interface_energy()
        dsasa = ia.get_interface_delta_sasa()
        nres = ia.get_num_interface_residues()
        dens = ddg / dsasa * 100.0 if dsasa else ""
        return ddg, dsasa, dens, cms(pose), nres
    finally:
        if tmp:
            try: os.unlink(tmp)
            except OSError: pass

def main():
    out = sys.argv[1]
    n_ours = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    w = open(out, "w")
    w.write("set,label,ddG,dSASA,ddG_per100,cms,nres\n"); w.flush()
    for d in sorted(glob.glob(DATA + "/*/")):
        cs = glob.glob(os.path.join(d, "afm3", "seed_*", "model.cif"))
        if not cs:
            continue
        try:
            ddg, dsasa, dens, c, nr = score(sorted(cs)[0])
            w.write("dataset,%s,%.2f,%.1f,%s,%s,%s\n" % (
                os.path.basename(d.rstrip("/")), ddg, dsasa,
                ("%.3f" % dens) if dens != "" else "", c, nr)); w.flush()
        except Exception as e:
            sys.stderr.write("ERR ds %s: %r\n" % (d, e))
    ours = glob.glob(RUN + "/chunk*/out/*_unrelaxed_rank_001*.pdb")
    random.seed(0); random.shuffle(ours); ours = ours[:n_ours]
    for f in ours:
        cand = os.path.basename(f).split("_unrelaxed")[0]
        try:
            ddg, dsasa, dens, c, nr = score(f)
            w.write("ours,%s,%.2f,%.1f,%s,%s,%s\n" % (
                cand, ddg, dsasa, ("%.3f" % dens) if dens != "" else "", c, nr)); w.flush()
        except Exception as e:
            sys.stderr.write("ERR ours %s: %r\n" % (cand, e))
    w.close()

if __name__ == "__main__":
    main()
