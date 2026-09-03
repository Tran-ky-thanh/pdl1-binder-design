#!/bin/bash
# Detailed progress view for the PD-L1 complex fold.
R=/home/thanh/protein_designs/pdl1_bench/fold_screen/complex_run
TOTAL=1050

bar() { # $1=done $2=total $3=width
  local d=$1 t=$2 w=${3:-28} fill pct
  [ "$t" -eq 0 ] && t=1
  fill=$(( d * w / t ))
  pct=$(( d * 100 / t ))
  printf '['
  printf '%0.s#' $(seq 1 $fill 2>/dev/null)
  printf '%0.s.' $(seq 1 $((w-fill)) 2>/dev/null)
  printf '] %3d%%' "$pct"
}

hhmm() { printf '%dh%02dm' $(( $1/3600 )) $(( ($1%3600)/60 )); }

# --- running? ---
if pgrep -f run_complex.sh >/dev/null || pgrep -f colabfold_batch >/dev/null; then
  echo "STATUS : RUNNING"
else
  if grep -q ALL_DONE "$R/run.log" 2>/dev/null; then echo "STATUS : FINISHED"; else echo "STATUS : NOT RUNNING (stopped/crashed)"; fi
fi

# --- elapsed since start ---
start=$(grep -m1 'COMPLEX_RUN START' "$R/run.log" 2>/dev/null | sed 's/.*START //; s/  host.*//')
if [ -n "$start" ]; then
  s0=$(date -d "$start" +%s 2>/dev/null); now=$(date +%s)
  echo "ELAPSED: $(hhmm $((now-s0)))  (since $start)"
fi

# --- overall + per-chunk ---
echo
gdone=0
for i in 0 1 2 3; do
  exp=$(ls "$R/chunk$i/in"/*.a3m 2>/dev/null | wc -l)
  don=$(ls "$R/chunk$i/out"/*.done.txt 2>/dev/null | wc -l)
  gdone=$((gdone+don))
  lb=$(sed -n "s/^${i},[^,]*,[^,]*,//p" "$R/manifest.csv" 2>/dev/null | sort -n | sed -n '1p;$p' | paste -sd- -)
  printf 'chunk%s  %3d/%-3d  %s  Lb %s\n' "$i" "$don" "$exp" "$(bar $don $exp 22)" "$lb"
done
echo "-----"
printf 'TOTAL   %4d/%-4d %s\n' "$gdone" "$TOTAL" "$(bar $gdone $TOTAL 28)"

# --- rate + ETA from done-marker timestamps ---
mapfile -t T < <(find "$R"/chunk*/out -name '*.done.txt' -printf '%T@\n' 2>/dev/null | sort -n)
n=${#T[@]}
if [ "$n" -ge 2 ]; then
  first=${T[0]%.*}; last=${T[$((n-1))]%.*}
  span=$(( last - first )); [ "$span" -lt 1 ] && span=1
  avg=$(( span / (n-1) ))                      # overall sec/complex
  k=10; [ "$n" -lt 11 ] && k=$((n-1))
  recent=$(( (last - ${T[$((n-1-k))]%.*}) / k ))   # recent sec/complex (last k)
  rem=$(( TOTAL - gdone ))
  eta=$(( rem * recent ))
  echo
  printf 'RATE   : %ds/complex (recent %d)   %d done, %d left\n' "$avg" "$recent" "$gdone" "$rem"
  printf 'ETA    : ~%s  (finish ~%s)\n' "$(hhmm $eta)" "$(date -d "@$(( $(date +%s) + eta ))" '+%a %H:%M')"
else
  echo; echo "RATE   : (need >=2 completed to estimate; warming up)"
fi

# --- currently folding ---
lastlog=$(ls -t "$R"/logs_chunk*.log 2>/dev/null | head -1)
if [ -n "$lastlog" ]; then
  cur=$(grep -aE 'Query [0-9]+/[0-9]+' "$lastlog" | tail -1 | sed 's/.*Query/Query/')
  rk=$(grep -aE 'rank_00.*pLDDT' "$lastlog" | tail -1 | sed 's/.*took//; s/^/last model /')
  echo
  echo "NOW    : ${cur:-<initializing / compiling>}   [$(basename "$lastlog" .log)]"
  [ -n "$rk" ] && echo "         $rk"
fi

# --- GPU ---
echo
echo -n "GPU    : "; nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu,temperature.gpu --format=csv,noheader 2>/dev/null
