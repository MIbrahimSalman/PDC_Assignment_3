#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
# Use >=10^8 paths even at p=1. Calibrated target is twelve seconds.
for p in 1 8; do
  t=$(awk -F '\t' -v p="$p" '$1=="mc_good" && $2==p {split($4,a," "); split(a[1],b,"="); print b[2]}' results/raw.tsv | sort -n | sed -n '2p')
  samples=$(awk -v t="$t" 'BEGIN {m=12*100000000/t; if(m<100000000)m=100000000; printf "%.0f", m}')
  for repetition in 1 2 3; do
    printf -v key 'mc_fixed_time12\t%s\t%s\t' "$p" "$repetition"
    if grep -Fq "$key" results/raw.tsv; then continue; fi
    echo "Running twelve-second Monte Carlo p=$p repetition=$repetition M=$samples"
    output=$(OMP_NUM_THREADS="$p" cs3_montecarlo/parallel "$samples" 42 good noverify)
    printf '%s%s\n' "$key" "$output" >> results/raw.tsv
  done
done
