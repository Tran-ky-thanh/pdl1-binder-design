#!/usr/bin/env bash
# =============================================================================
# PD-L1 de novo binder pipeline - single entry point / orchestrator.
#
#   ./run_pipeline.sh list            # show all stages
#   ./run_pipeline.sh check           # prereq / environment check for every stage
#   ./run_pipeline.sh <stage>         # run (or document) one stage, e.g. 07_rank
#   ./run_pipeline.sh all             # run local stages, document the GPU/cloud ones
#
# Heavy stages (RFdiffusion, MPNN, ColabFold, cloud co-folders) need external
# tools / GPU / API keys and are only *documented + prereq-checked* unless their
# tools are present and RUN_HEAVY=1. The local stages (score, rank, report) run
# for real. Full write-up, per-stage tips and env details: PIPELINE.md
#
# Two environments are used - the script prints the exact switch command each
# time it needs one, so you never have to guess:
#   * ROOT  env : RFdiffusion + localColabFold, installed under /root/protein
#                 -> run as root, conda at /root/miniforge3   (stages 00,02,04)
#   * USER venv : scoring / ranking / analysis, $PROJECT/.venv (stages 05,07,08)
#                 numpy is toggled 1.26.4 (DockQ) <-> 2.5.2 (ESM) - see ensure_numpy
# =============================================================================
set -uo pipefail

# ---------- config (override via environment) ----------
PROJECT="${PROJECT:-/home/thanh/protein_designs/pdl1_bench}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$REPO/scripts"
VENV="${VENV:-$PROJECT/.venv}"
PY="$VENV/bin/python"
PIP="$VENV/bin/pip"
SCORING="$PROJECT/fold_screen/complex_run/scoring"
COMPLEX="$PROJECT/fold_screen/complex_run"
RUN_HEAVY="${RUN_HEAVY:-0}"
DRYRUN="${DRYRUN:-0}"

# external tool roots (run as root)
CONDA_ROOT="/root/miniforge3"
RFDIFF="/root/protein/RFdiffusion"
LCF="/root/protein/localcolabfold"

# ColabFold time-savers (see PIPELINE.md, stage 04)
export JAX_COMPILATION_CACHE_DIR="${JAX_COMPILATION_CACHE_DIR:-$HOME/.cache/colabfold_jax}"
TARGET_MSA="$PROJECT/fold_screen/target_msa"      # PD-L1 MSA computed ONCE, reused for every binder

# cloud co-folder credentials - MUST come from the environment, never hardcoded
: "${BOLTZ_API_KEY:=}"        # api.boltz.bio     (live key sk_bc_ws_live_...)
: "${BIOHUB_TOKEN:=}"         # ESMFold2-Full via Biohub Forge
# Protenix (JapanFold) needs no key.

# ---------- pretty ----------
c_hdr(){ printf '\n\033[1;36m== %s ==\033[0m\n' "$*"; }
c_env(){ printf '   \033[1;33mENV \033[0m %s\n' "$*"; }
c_run(){ printf '   \033[1;32mRUN \033[0m %s\n' "$*"; }
note(){ printf '   \033[2m%s\033[0m\n' "$*"; }
tip(){  printf '   \033[35mTIP \033[0m %s\n' "$*"; }
ok(){   printf '   \033[1;32mOK  \033[0m %s\n' "$*"; }
skip(){ printf '   \033[1;31mSKIP\033[0m %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }
run(){ c_run "$*"; [ "$DRYRUN" = 1 ] && return 0; eval "$*"; }

# ---------- env switches (printed AND applied) ----------
print_root_env(){
  c_env "enter root + activate localColabFold/RFdiffusion conda:"
  note   "sudo -i"
  note   "source $CONDA_ROOT/etc/profile.d/conda.sh"
  note   "export PATH=$LCF/colabfold-conda/bin:\$PATH"
  note   "export JAX_COMPILATION_CACHE_DIR=$JAX_COMPILATION_CACHE_DIR ; unset DISPLAY ; export MPLBACKEND=Agg"
}
ensure_numpy(){   # $1 = required version; toggles the venv's numpy if needed
  local want="$1" cur
  cur="$($PY -c 'import numpy;print(numpy.__version__)' 2>/dev/null)"
  if [ "$cur" != "$want" ]; then
    c_env "numpy $cur -> $want  (DockQ needs <2, ESM SDK needs 2.x)"
    note "$PIP install -q numpy==$want"
    [ "$DRYRUN" = 1 ] || $PIP install -q "numpy==$want"
  else
    c_env "numpy already $want"
  fi
}

find_complex_pdb(){  # $1 = cand -> echo path of AF2-multimer complex, or empty
  local c="$1" g
  for i in 0 1 2 3; do
    g=$(ls "$COMPLEX"/chunk$i/out/${c}_unrelaxed_rank_001_*.pdb 2>/dev/null | head -1)
    [ -n "$g" ] && { echo "$g"; return; }
  done
}

# =============================================================================
# STAGES
# =============================================================================
stage_00_rfdiff(){
  c_hdr "00  RFdiffusion  ->  52 backbones"
  tip "Seed from real epitope: crystal complexes 5C3T/4ZQK/4Z18 + de novo scaffolds,"
  tip "hotspots Y56/R113/M115/Y123 fixed. Anchoring on a known epitope keeps the whole funnel on a bindable face."
  print_root_env
  c_run "python $RFDIFF/scripts/run_inference.py \\"
  note  "  inference.output_prefix=$PROJECT/rfdiff_out2/4ZQK_binder \\"
  note  "  inference.input_pdb=<PD-L1_target.pdb> inference.num_designs=... \\"
  note  "  'contigmap.contigs=[A1-115/0 60-90]' 'ppi.hotspot_res=[A56,A113,A115,A123]' \\"
  note  "  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0   # checkpoint: Complex_base_ckpt.pt"
  [ -d "$RFDIFF" ] && ok "RFdiffusion present: $RFDIFF" || skip "RFdiffusion not found (needs root + GPU) - documented only"
}

stage_01_mpnn(){
  c_hdr "01  ProteinMPNN + SolubleMPNN  ->  4,905 sequences"
  tip "Run BOTH sequence designers at sampling temp 0.2. SolubleMPNN folds better as a monomer"
  tip "(median pLDDT 93.0 vs 90.3) - it biases to soluble/stable seqs; not yet proof of better binding."
  c_env "ProteinMPNN's own env (its repo venv/conda); CPU is fine."
  c_run "python ProteinMPNN/protein_mpnn_run.py --pdb_path_chains B --sampling_temp 0.2 --num_seq_per_target N ..."
  c_run "python ProteinMPNN/protein_mpnn_run.py --use_soluble_model  --sampling_temp 0.2 ...   # SolubleMPNN"
  note "output parsed + per-backbone capped -> $PROJECT/mpnn_out/filtered.csv (2,481 candidate pool)"
  [ -f "$PROJECT/mpnn_out/filtered.csv" ] && ok "MPNN output present" || skip "run MPNN first"
}

stage_02_monomer(){
  c_hdr "02  localColabFold monomer foldability  ->  2,162 pass (pLDDT >= 70)"
  tip "Monomer foldability barely discriminates (pLDDT/pTM r 0.88-0.91, MPNN seqs fold well by construction);"
  tip "keep the gate LOOSE at pLDDT>=70 - over-filtering here throws away good binders for no signal."
  print_root_env
  c_run "colabfold_batch --num-recycle 1 --num-models 1 <seqs_fasta_dir> <out_dir>"
  tip "JAX compile cache ($JAX_COMPILATION_CACHE_DIR) is reused across runs - first fold compiles, rest are fast."
  [ -f "$PROJECT/fold_screen/gate_survivors.csv" ] && ok "monomer gate results present" || skip "run monomer fold first"
}

stage_03_diversity(){
  c_hdr "03  redundancy + diversity down-sample  ->  1,050"
  tip "Cheap filters first - but CHECK they bite: 90% dedup removed 0 here (diverse backbones, low MPNN temp)."
  tip "The real reducer is 70% diversity down-sampling, kept AFTER monomer folding so the cluster representative"
  tip "is the best-FOLDING one (needs pLDDT/pTM), not a blind pick by MPNN score."
  c_env "USER venv (mmseqs2 / cd-hit for clustering if available)"
  note "dedup 90% id -> 0 removed ; cluster 70% id, pick best-pLDDT per cluster -> 1,050"
  [ -f "$COMPLEX/manifest.csv" ] && ok "1,050 selection present ($COMPLEX/manifest.csv)" || skip "run diversity step first"
}

stage_04_complex(){
  c_hdr "04  localColabFold complex co-fold  ->  1,050 predicted complexes"
  tip "Reuse the target MSA: PD-L1's MSA is computed ONCE ($TARGET_MSA) and paired into every binder's .a3m,"
  tip "so the expensive MSA step runs once, not 1,050 times."
  tip "Batch it: split into 4 chunks (chunk0..3); the driver auto-resumes via .done.txt and restarts ColabFold on SIGKILL."
  tip "WSL: poll the log SPARINGLY (20-30 min) - aggressive wsl.exe polling that times out can restart the VM and wipe /tmp."
  print_root_env
  c_run "bash $SCRIPTS/run_complex.sh    # recycle 1, 1 model, reused MSA, 4-chunk resumable loop"
  if [ "$RUN_HEAVY" = 1 ] && have colabfold_batch; then run "bash $SCRIPTS/run_complex.sh"; else skip "GPU/root stage - set RUN_HEAVY=1 as root to execute"; fi
}

stage_05_score(){
  c_hdr "05  interface scoring (local, CPU)   [RUNNABLE]"
  tip "ipSAE is min-of-both-directions: min of the A->B and B->A asym rows for ONE structure (not min across models);"
  tip "it is pure-PAE, so it is invariant to relaxation - safe to compute once on the raw fold."
  tip "Relax the BACKBONE before ddG: side-chain-only relax on recycle-1 structures invents clash blow-ups;"
  tip "use a coordinate-constrained FastRelax (backbone free). ipSAE/sc_DockQ barely move; energy terms become usable."
  ensure_numpy 1.26.4
  local cand="${1:-cand00169}" mp; mp="$(find_complex_pdb "$cand")"
  if [ -z "$mp" ]; then skip "no complex PDB for $cand (run stage 04 first)"; return; fi
  local scores; scores="$(ls "$(dirname "$mp")"/${cand}_scores_rank_001_*.json 2>/dev/null | head -1)"
  c_run "$PY $SCRIPTS/ipsae.py <pae.json> <complex.pdb> 10 10       # -> ipSAE_min"
  [ -n "$scores" ] && run "$PY $SCRIPTS/ipsae.py '$scores' '$mp' 10 10 >/dev/null 2>&1 && echo '   ipSAE ran on $cand'"
  c_run "$PY $SCRIPTS/score_iface.py $cand <complex.pdb>           # -> sc_DockQ,BSA,contacts"
  run "$PY $SCRIPTS/score_iface.py '$cand' '$mp'"
  c_run "$PY $SCRIPTS/bench_ddg.py ...   # PyRosetta ddG/CMS after CONSTRAINED FULL relax (needs pyrosetta)"
}

stage_06_cofold(){
  c_hdr "06  cloud multi-predictor consensus  ->  43 pose-PASS"
  tip "Three opinions beat two: Boltz-2 + ESMFold2-Full agreed on 46; adding Protenix rejected 3 more"
  tip "(high Boltz/ESM sc_DockQ but Protenix <0.23). A cheap third co-folder is the highest-value filter here."
  tip "Free co-folders: ESMFold2-Full (Biohub Forge, PAE-capable -> ipSAE works); Protenix (JapanFold, NO PAE -> iptm+sc_DockQ only);"
  tip "Boltz-2 test keys return a MOCK fixture - only live keys fold real sequences."
  c_env "For ESMFold2 SDK you need numpy 2.x:"; ensure_numpy_hint
  note "BOLTZ_API_KEY $( [ -n "$BOLTZ_API_KEY" ] && echo SET || echo 'MISSING (export it)') | BIOHUB_TOKEN $( [ -n "$BIOHUB_TOKEN" ] && echo SET || echo 'MISSING (export it)')"
  c_run "$PY $SCRIPTS/cofold/boltz2_run.py    --key \$BOLTZ_API_KEY  --in cands.fasta --out scoring/boltz_ipsae_scdockq.csv"
  c_run "$PY $SCRIPTS/cofold/esmfold2_full.py --token \$BIOHUB_TOKEN --in cands.fasta --out scoring/esm_full_ipsae_scdockq.csv"
  c_run "$PY $SCRIPTS/cofold/protenix_run.py                        --in cands.fasta --out scoring/protenix_results.csv"
  note "(these API clients lived in /tmp and were lost to a VM restart - re-create from the outlines in PIPELINE.md)"
  skip "cloud stage - needs API keys; documented in PIPELINE.md"
}
ensure_numpy_hint(){ note "$PIP install -q numpy==2.5.2   # for the ESM SDK; switch back to 1.26.4 for DockQ scoring"; }

stage_07_rank(){
  c_hdr "07  rank designs  [RUNNABLE]"
  tip "Composite z = 4*z(ipSAE) + 1*z(sc_DockQ), the dataset's 4:1 weighting, over the pose-passing designs."
  ensure_numpy 1.26.4
  c_run "$PY $SCRIPTS/analysis/rank_designs.py --scoring $SCORING --verify"
  run "$PY $SCRIPTS/analysis/rank_designs.py --scoring '$SCORING' --verify"
}

stage_08_report(){
  c_hdr "08  build interactive report"
  note "superpose top-30 on PD-L1, slim the data, inject into the template:"
  c_run "$PY $SCRIPTS/analysis/novelty_and_extract.py   # novelty vs dataset + extract complexes -> report_data.json"
  c_run "$PY $SCRIPTS/analysis/slim_report_data.py      # -> report_slim.json (top-30 structures)"
  c_run "$PY $SCRIPTS/analysis/align_top30.py           # superpose all top-30 onto the PD-L1 frame"
  c_run "$PY $SCRIPTS/analysis/build_report.py          # inject data into report.template.html -> report.html"
  note "prebuilt report is committed at report/index.html (open in a browser; 3Dmol.js from CDN)."
  note "these analysis scripts currently use absolute scratch paths - see PIPELINE.md before re-running."
  [ -f "$REPO/report/index.html" ] && ok "report present: report/index.html" || skip "report not built yet"
}

STAGES=(00_rfdiff 01_mpnn 02_monomer 03_diversity 04_complex 05_score 06_cofold 07_rank 08_report)

do_check(){
  c_hdr "environment check"
  printf '   PROJECT = %s\n   REPO    = %s\n   VENV    = %s\n' "$PROJECT" "$REPO" "$VENV"
  [ -x "$PY" ] && ok "venv python: $($PY --version 2>&1)" || skip "venv python missing at $PY"
  $PY -c 'import Bio' 2>/dev/null && ok "biopython" || skip "biopython missing"
  $PY -c 'import DockQ' 2>/dev/null && ok "DockQ" || skip "DockQ missing (pip install DockQ==2.1.3, numpy<2)"
  $PY -c 'import pyrosetta' 2>/dev/null && ok "pyrosetta" || skip "pyrosetta missing (ddG stage only)"
  [ -d "$RFDIFF" ] && ok "RFdiffusion $RFDIFF" || skip "RFdiffusion (root/GPU) not found"
  [ -d "$LCF" ] && ok "localColabFold $LCF" || skip "localColabFold (root/GPU) not found"
  [ -d "$TARGET_MSA" ] && ok "reusable target MSA $TARGET_MSA" || note "target MSA cache absent"
  [ -d "$JAX_COMPILATION_CACHE_DIR" ] && ok "JAX compile cache $JAX_COMPILATION_CACHE_DIR" || note "JAX cache absent (built on first fold)"
  [ -d "$SCORING" ] && ok "scoring dir $SCORING" || skip "scoring dir absent"
  note "keys: BOLTZ_API_KEY $( [ -n "$BOLTZ_API_KEY" ] && echo SET || echo unset ) | BIOHUB_TOKEN $( [ -n "$BIOHUB_TOKEN" ] && echo SET || echo unset )"
}

main(){
  local a="${1:-list}"
  case "$a" in
    list)  c_hdr "stages"; for s in "${STAGES[@]}"; do echo "   $s"; done
           note "run one: ./run_pipeline.sh 07_rank   |  all: ./run_pipeline.sh all  |  ./run_pipeline.sh check";;
    check) do_check;;
    all)   for s in "${STAGES[@]}"; do "stage_$s"; done;;
    *)     if printf '%s\n' "${STAGES[@]}" | grep -qx "$a"; then shift; "stage_$a" "$@";
           else echo "unknown stage '$a'"; echo "try: ./run_pipeline.sh list"; exit 2; fi;;
  esac
}
main "$@"
