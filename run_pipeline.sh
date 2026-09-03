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
PROJECT="${PROJECT:-$HOME/protein_designs/pdl1_bench}"
# Python scripts read these; export so subprocesses (scoring, cofold) inherit them.
export PDL1_PROJECT="${PDL1_PROJECT:-$PROJECT}"
export PDL1_OUT="${PDL1_OUT:-$PROJECT/build}"
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
RFDIFF="${RFDIFF:-/root/protein/RFdiffusion}"
LCF="${LCF:-/root/protein/localcolabfold}"

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
  note "inputs committed at inputs/rfdiffusion/ (pdl1_target/5C3T/4ZQK/4Z18). 4 arms, exact contigs+hotspots:"
  c_run "bash $SCRIPTS/00_rfdiffusion.sh          # de novo A2-110/0 70-80 + crystal A18-132/0 65-85"
  note "hotspots: crystal A56,A113,A115,A123  |  de-novo A39,A96,A98,A106  (same epitope, +17 offset); noise=0"
  [ -d "$RFDIFF" ] && ok "RFdiffusion present: $RFDIFF" || skip "RFdiffusion not found (needs root + GPU) - documented only"
  [ -d "$REPO/inputs/rfdiffusion" ] && ok "target inputs present: inputs/rfdiffusion/" || note "inputs/rfdiffusion missing"
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
  note "Score ALL 1,050 co-folds (ipSAE + sc_DockQ + PyRosetta ΔΔG/CMS) -> the ranked_1050 table that 05b consumes:"
  c_run "$PY $SCRIPTS/score_complexes.py --out \$PDL1_OUT/iface_1050.csv   # needs PyRosetta + the 1,050 complexes; parallelise with --slice K N"
  note "bench_ddg.py below is only a QC benchmark (dataset vs ours on ~50); the per-design smoke test follows:"
  local cand="${1:-cand00169}" mp; mp="$(find_complex_pdb "$cand")"
  if [ -z "$mp" ]; then skip "no complex PDB for $cand (run stage 04 first)"; return; fi
  local scores; scores="$(ls "$(dirname "$mp")"/${cand}_scores_rank_001_*.json 2>/dev/null | head -1)"
  c_run "$PY $SCRIPTS/ipsae.py <pae.json> <complex.pdb> 10 10       # -> ipSAE_min"
  [ -n "$scores" ] && run "$PY $SCRIPTS/ipsae.py '$scores' '$mp' 10 10 >/dev/null 2>&1 && echo '   ipSAE ran on $cand'"
  c_run "$PY $SCRIPTS/score_iface.py $cand <complex.pdb>           # -> sc_DockQ,BSA,contacts"
  run "$PY $SCRIPTS/score_iface.py '$cand' '$mp'"
  c_run "$PY $SCRIPTS/bench_ddg.py ...   # PyRosetta ddG/CMS after CONSTRAINED FULL relax (needs pyrosetta)"
  note "ddG/CMS are computed for all 1,050 (ranked_1050.csv) and feed the Stage-2 physics gate (see 05b_filter)."
}

stage_05b_filter(){
  c_hdr "05b  two-stage screen: ipSAE+sc_DockQ -> 242 ; PyRosetta physics -> 157   [RUNNABLE]"
  tip "Stage-2 is a PHYSICS gate, not another pose gate: all 242 Stage-1 designs already pass sc_DockQ>=0.23."
  tip "PyRosetta ΔΔG<=-40 REU & CMS>=360 Å² removes 85 weak/repulsive interfaces (242 -> 157); then curation -> 50."
  ensure_numpy 1.26.4
  c_run "$PY $SCRIPTS/analysis/filter_stages.py --ranked $REPO/results/ranked_1050.csv --stage1 $REPO/results/stage1_pass.txt --verify"
  run "$PY $SCRIPTS/analysis/filter_stages.py --ranked '$REPO/results/ranked_1050.csv' --stage1 '$REPO/results/stage1_pass.txt' --verify"
}

stage_06_cofold(){
  c_hdr "06  cloud multi-predictor consensus  ->  43 pose-PASS"
  tip "Three opinions beat two: Boltz-2 + ESMFold2-Full agreed on 46; adding Protenix rejected 3 more"
  tip "(high Boltz/ESM sc_DockQ but Protenix <0.23). A cheap third co-folder is the highest-value filter here."
  tip "Free co-folders: ESMFold2-Full (Biohub Forge, PAE-capable -> ipSAE works); Protenix (JapanFold, NO PAE -> iptm+sc_DockQ only);"
  tip "Boltz-2 test keys return a MOCK fixture - only live keys fold real sequences."
  c_env "For ESMFold2 SDK you need numpy 2.x:"; ensure_numpy_hint
  note "BOLTZ_API_KEY $( [ -n "$BOLTZ_API_KEY" ] && echo SET || echo 'MISSING (export it)') | BIOHUB_TOKEN $( [ -n "$BIOHUB_TOKEN" ] && echo SET || echo 'MISSING (export it)')"
  note "All three co-folders are RUNNABLE and tested. Input: results/final_selection.csv (cand,seq)."
  note "-- Boltz-2 (numpy 1.26.4; needs a LIVE key) --"
  c_run "$PY $SCRIPTS/cofold/boltz2_run.py --in $REPO/results/final_selection.csv --out $SCORING/boltz_ipsae_scdockq.csv --scoring $SCORING --scripts $SCRIPTS --key \$BOLTZ_API_KEY"
  note "-- ESMFold2-Full (TWO passes; the numpy majors cannot coexist) --"
  c_env "pip install -q numpy==2.5.2   # FOLD pass (esm SDK)"
  c_run "$PY $SCRIPTS/cofold/esmfold2_full.py fold  --in $REPO/results/final_selection.csv --outdir /tmp/esm_full_struct --results $SCORING/esmfold2_full_results.csv --token \$BIOHUB_TOKEN"
  c_env "pip install -q numpy==1.26.4  # SCORE pass (DockQ)"
  c_run "$PY $SCRIPTS/cofold/esmfold2_full.py score --outdir /tmp/esm_full_struct --out $SCORING/esm_full_ipsae_scdockq.csv --scoring $SCORING --scripts $SCRIPTS"
  note "-- Protenix (JapanFold; no key; batched for the rate limit) --"
  c_run "$PY $SCRIPTS/cofold/protenix_run.py --in $REPO/results/final_selection.csv --out $SCORING/protenix_results.csv --scoring $SCORING --scripts $SCRIPTS --batch 6 --space 3.5 --cooldown 150"
  if [ "$RUN_HEAVY" = 1 ]; then ensure_numpy 1.26.4; run "$PY $SCRIPTS/cofold/protenix_run.py --in $REPO/results/final_selection.csv --out $SCORING/protenix_results.csv --scoring $SCORING --scripts $SCRIPTS --batch 6 --space 3.5 --cooldown 150"; else skip "set RUN_HEAVY=1 to submit the Protenix batch; run Boltz/ESM manually (see commands above)"; fi
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
  c_run "$PY $SCRIPTS/analysis/novelty_and_extract.py   # novelty (sequence) vs dataset + extract complexes -> report_data.json"
  c_run "$PY $SCRIPTS/analysis/slim_report_data.py      # -> report_slim.json (top-30 structures)"
  c_run "$PY $SCRIPTS/analysis/align_top30.py           # superpose all top-30 onto the PD-L1 frame"
  note "structural novelty (TM-align) needs the full Anthropic PD-L1 reference set and a TMalign binary:"
  c_run "$PY $SCRIPTS/analysis/fetch_pdl1_designs.py    # download 'designed' backbones for all 90 PD-L1 rows (78/90 have one)"
  note "TMalign is a standalone C++ binary (not pip-installable): g++ -O2 -o ~/bin/TMalign TMalign.cpp"
  note "  source: https://zhanggroup.org/TM-align/TMalign.cpp (Zhang Lab); Foldseek was tried first but its"
  note "  150 MB release binary was infeasible over a throttled connection, so TM-align (~200 KB source) was used instead."
  c_run "$PY $SCRIPTS/analysis/structural_novelty.py --tmalign ~/bin/TMalign --top 30  # -> results/tm_vs_anthropic.csv"
  c_run "$PY $SCRIPTS/analysis/make_valfig.py           # dataset-consensus validation figure -> colabfold_consensus.png"
  c_run "$PY $SCRIPTS/analysis/build_extras.py          # embed figures + novelty/fold-pose 3Dmol overlays into report_slim.json"
  c_run "$PY $SCRIPTS/analysis/build_report.py          # inject data into report.template.html -> report.html"
  c_run "cp \$PDL1_OUT/report.html $REPO/report/index.html"
  note "prebuilt report is committed at report/index.html (open in a browser; 3Dmol.js from CDN)."
  note "set PDL1_PROJECT and PDL1_OUT first (see PIPELINE.md); generated files land in \$PDL1_OUT."
  [ -f "$REPO/report/index.html" ] && ok "report present: report/index.html" || skip "report not built yet"
}

STAGES=(00_rfdiff 01_mpnn 02_monomer 03_diversity 04_complex 05_score 05b_filter 06_cofold 07_rank 08_report)

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
