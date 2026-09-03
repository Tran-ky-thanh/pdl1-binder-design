#!/usr/bin/env python3
"""OPTIONAL stronger layer: PyRosetta interface energetics (ddG / CMS / SC).
Baker-lab de-novo minibinder filter metrics (Cao et al. 2022, Nature).
Requires a one-time PyRosetta install (academic license, free):
    pip install pyrosetta-installer
    python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"

Per candidate (AF2 complex, chain A=target, B=binder):
  ddG        : dG_separated from InterfaceAnalyzer (binding energy, REU; lower=better)
  dSASA      : interface buried SASA (A^2)
  ddG_dSASA  : ddG per 100 A^2 (energy density; the key minibinder discriminator)
  sc         : shape complementarity (0-1)
  cms        : contact molecular surface (if ContactMolecularSurfaceFilter available)
  nres_int   : interface residue count
FastRelax is applied (bb+sc, 1 repeat) before analysis for stable ddG; set
RELAX=0 to score raw model (faster, noisier).

Usage: python3 rosetta_iface.py CAND MODEL_PDB
       -> cand,ddG,dSASA,ddG_dSASA,sc,cms,nres_int
UNTESTED at authoring time (no PyRosetta) -- validate on 2-3 complexes first.
"""
import os, sys
import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

pyrosetta.init("-mute all -use_input_sc -ex1 -ex2aro")

def relax(pose):
    from pyrosetta.rosetta.protocols.relax import FastRelax
    sf = pyrosetta.get_fa_scorefxn()
    fr = FastRelax(sf, 1)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(True); mm.set_chi(True)
    fr.set_movemap(mm)
    fr.apply(pose)

def _chain_sel(ch):
    return pyrosetta.rosetta.core.select.residue_selector.ChainSelector(ch)

def cms(pose):
    try:
        F = pyrosetta.rosetta.protocols.simple_filters.ContactMolecularSurfaceFilter
        f = F(); f.selector1(_chain_sel("A")); f.selector2(_chain_sel("B"))
        return round(f.compute(pose), 1)
    except Exception:
        return ""

def main():
    cand, pdb = sys.argv[1], sys.argv[2]
    pose = pose_from_pdb(pdb)
    if os.environ.get("RELAX", "1") == "1":
        relax(pose)
    ia = InterfaceAnalyzerMover("A_B")
    ia.set_compute_packstat(True)
    ia.set_pack_separated(True)
    ia.apply(pose)
    data = ia.get_all_data()
    ddg = ia.get_separated_interface_energy()
    dsasa = ia.get_interface_delta_sasa()
    sc = getattr(data, "sc_value", "") if data else ""
    nres = ia.get_num_interface_residues()
    dens = (ddg / dsasa * 100.0) if dsasa else ""
    print("%s,%.2f,%.1f,%s,%s,%s,%s" % (
        cand, ddg, dsasa,
        ("%.3f" % dens) if dens != "" else "",
        ("%.3f" % sc) if sc not in ("", None) else "",
        cms(pose), nres))

if __name__ == "__main__":
    main()
