#!/usr/bin/env bash
# =============================================================================
# Stage 00 - RFdiffusion binder-backbone generation for PD-L1.  (ROOT env, GPU)
#
#   sudo -i
#   source /root/miniforge3/etc/profile.d/conda.sh
#   conda activate <your RFdiffusion env>     # SE3nv / rfdiffusion
#   bash 00_rfdiffusion.sh
#
# Reproduces the 4 generation arms used in this project. Inputs are the trimmed
# PD-L1 target structures in ../inputs/rfdiffusion/. Backbones land in $OUT, then
# feed stage 01 (ProteinMPNN/SolubleMPNN). Checkpoint auto-selected: Complex_base_ckpt.pt
#
# Exact settings recovered from the run .trb configs:
#   arm           input_pdb        contig             binder len  hotspot_res (author numbering)
#   pdl1_binder   pdl1_target.pdb  A2-110/0  70-80    70-80       A39,A96,A98,A106
#   5C3T_binder   pdl1_5C3T.pdb    A18-132/0 65-85    65-85       A56,A113,A115,A123
#   4ZQK_binder   pdl1_4ZQK.pdb    A18-132/0 65-85    65-85       A56,A113,A115,A123
#   4Z18_binder   pdl1_4Z18.pdb    A18-132/0 65-85    65-85       A56,A113,A115,A123
# (pdl1_target is numbered 2-110; the crystal templates 18-132 -> the same epitope,
#  hotspots differ only by the +17 numbering offset. Zero noise = deterministic-ish backbones.)
# =============================================================================
set -euo pipefail

RFDIFF="${RFDIFF:-/root/protein/RFdiffusion}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IN="${IN:-$HERE/../inputs/rfdiffusion}"
OUT="${OUT:-$HOME/rfdiff_out}"
NUM="${NUM:-8}"          # designs per invocation; we looped this to build ~50 backbones/template
mkdir -p "$OUT"

gen(){  # $1 name  $2 input_pdb  $3 contig  $4 hotspots
  echo ">> $1  ($2  contig=[$3]  hotspots=[$4])"
  python "$RFDIFF/scripts/run_inference.py" \
    inference.input_pdb="$2" \
    inference.output_prefix="$OUT/$1" \
    inference.num_designs="$NUM" \
    "contigmap.contigs=[$3]" \
    "ppi.hotspot_res=[$4]" \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0
}

gen pdl1_binder "$IN/pdl1_target.pdb" "A2-110/0 70-80"  "A39,A96,A98,A106"
gen 5C3T_binder "$IN/pdl1_5C3T.pdb"   "A18-132/0 65-85" "A56,A113,A115,A123"
gen 4ZQK_binder "$IN/pdl1_4ZQK.pdb"   "A18-132/0 65-85" "A56,A113,A115,A123"
gen 4Z18_binder "$IN/pdl1_4Z18.pdb"   "A18-132/0 65-85" "A56,A113,A115,A123"

echo ">> done. Backbones in $OUT ; next: stage 01 (ProteinMPNN + SolubleMPNN)."
