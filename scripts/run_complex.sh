#!/bin/bash
# PD-L1 complex fold — 1050 binders vs PD-L1, AF2-multimer, recycle 1, reused MSA.
# RESILIENT: colabfold has been observed to die (SIGKILL/137) after hours on WSL2.
# This loop always attacks the lowest incomplete chunk and auto-restarts on death,
# grinding until all chunks complete. .done.txt persistence => every restart resumes.
source /root/miniforge3/etc/profile.d/conda.sh
export PATH=/root/protein/localcolabfold/colabfold-conda/bin:$PATH
export MPLBACKEND=Agg
unset DISPLAY
export JAX_COMPILATION_CACHE_DIR=$HOME/.cache/colabfold_jax
mkdir -p "$JAX_COMPILATION_CACHE_DIR"

ROOT=/home/thanh/protein_designs/pdl1_bench/fold_screen/complex_run
cd "$ROOT"

exec 9>"$ROOT/.run.lock"
if ! flock -n 9; then
  echo "ALREADY_RUNNING skip $(date)" >> run.log
  exit 0
fi

TOTAL=1050
count_done(){ find "$ROOT"/chunk*/out -name '*.done.txt' 2>/dev/null | wc -l; }
echo "COMPLEX_RUN START $(date)  host=$(hostname)  resume_from=$(count_done)/$TOTAL" >> run.log

stall=0; pass=0
while :; do
  pass=$((pass+1))
  before=$(count_done)

  target=-1
  for i in 0 1 2 3; do
    exp=$(ls "$ROOT/chunk$i/in"/*.a3m 2>/dev/null | wc -l)
    don=$(ls "$ROOT/chunk$i/out"/*.done.txt 2>/dev/null | wc -l)
    if [ "$don" -lt "$exp" ]; then target=$i; break; fi
  done

  if [ "$target" -lt 0 ]; then
    echo "COMPLEX_RUN ALL_DONE $(count_done)/$TOTAL $(date)" >> run.log
    break
  fi

  i=$target
  ind="$ROOT/chunk$i/in"; outd="$ROOT/chunk$i/out"; mkdir -p "$outd"
  exp=$(ls "$ind"/*.a3m 2>/dev/null | wc -l)
  don=$(ls "$outd"/*.done.txt 2>/dev/null | wc -l)
  echo "PASS $pass START chunk$i ($don/$exp) total=$(count_done)/$TOTAL $(date +%H:%M:%S)" >> run.log
  t0=$(date +%s)
  colabfold_batch --num-recycle 1 --num-models 1 \
    --max-seq 256 --max-extra-seq 512 --sort-queries-by length \
    "$ind" "$outd" >> "$ROOT/logs_chunk$i.log" 2>&1
  rc=$?
  t1=$(date +%s)
  don2=$(ls "$outd"/*.done.txt 2>/dev/null | wc -l)
  echo "PASS $pass END chunk$i $((t1-t0))s exit $rc done=$don2/$exp total=$(count_done)/$TOTAL $(date +%H:%M:%S)" >> run.log

  after=$(count_done)
  if [ "$after" -le "$before" ]; then
    stall=$((stall+1))
    echo "  NO_PROGRESS stall=$stall rc=$rc $(date +%H:%M:%S)" >> run.log
    if [ "$stall" -ge 6 ]; then
      echo "COMPLEX_RUN ABORT stalled 6x no progress $(date)" >> run.log
      break
    fi
    sleep 30
  else
    stall=0
    sleep 3
  fi
done