#!/bin/bash
# CPU orthogonal screen driver (sc_DockQ + interface geometry). NO GPU.
# Resumable (per-cand .row files) + parallel. Safe to run while GPU folds.
# Usage: bash run_scoring.sh [NPROC]   (default 4)
set -u
RUN=/home/thanh/protein_designs/pdl1_bench/fold_screen/complex_run
SC=$RUN/scoring; ROWS=$SC/rows; mkdir -p "$ROWS"
PY=/home/thanh/protein_designs/pdl1_bench/.venv/bin/python
SCRIPT=$SC/score_iface.py
NPROC=${1:-4}
cd "$SC"

: > worklist.txt
for f in $RUN/chunk*/out/*_unrelaxed_rank_001*.pdb; do
  [ -e "$f" ] || continue
  cand=$(basename "$f" | sed 's/_unrelaxed.*//')
  [ -s "$ROWS/$cand.row" ] && continue
  printf '%s\t%s\n' "$cand" "$f" >> worklist.txt
done
n=$(wc -l < worklist.txt)
echo "to score: $n  (parallel=$NPROC)  $(date +%H:%M:%S)"
[ "$n" -eq 0 ] || xargs -P "$NPROC" -a worklist.txt -L1 bash -c '
  "'"$PY"'" "'"$SCRIPT"'" "$0" "$1" > "'"$ROWS"'/$0.row" 2>"'"$ROWS"'/$0.err" \
    || echo "$0,ERR" > "'"$ROWS"'/$0.row"
'
echo "cand,sc_dockq,irms,Lrms,fnat,bsa,n_int_res,n_contact" > iface_scores.csv
cat "$ROWS"/*.row >> iface_scores.csv 2>/dev/null
echo "aggregated $(( $(wc -l < iface_scores.csv) - 1 )) rows -> $SC/iface_scores.csv  $(date +%H:%M:%S)"
