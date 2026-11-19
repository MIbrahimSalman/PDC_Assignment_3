#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
for p in 1 2 4; do
  export OMP_NUM_THREADS="$p"
  for samples in 2 65537 1000003; do
    cs3_montecarlo/parallel "$samples" 42 good verify
  done
  for edges in 0 3 10003; do
    for mode in pull push; do
      for specification in 'static 0' 'static 3' 'dynamic 1' 'dynamic 17' 'guided 4'; do
        read -r schedule chunk <<< "$specification"
        cs4_pagerank_openmp/parallel 128 "$edges" 42 "$mode" "$schedule" "$chunk" verify
      done
    done
  done
  for n in 1 3 10003; do
    for distribution in uniform skew; do
      for oversampling in 1 8; do
        mpirun -np "$p" cs5_samplesort_mpi/parallel "$n" "$distribution" "$oversampling" 42 verify
      done
    done
  done
done
echo 'Remaining cases passed numerical/exact serial comparisons.'
