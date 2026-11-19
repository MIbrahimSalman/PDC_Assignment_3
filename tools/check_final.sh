#!/usr/bin/env bash
# Additional independent invariants and non-square six-rank layouts.
# Run after benchmarks, not concurrently with measurements.
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=6 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
cs1_heat_mpi/serial 3 1 1 | awk '
  {for(i=1;i<=NF;i++){split($i,a,"=");v[a[1]]=a[2]}}
  END {if(v["checksum"]!=320 || v["max_change"]!=20) exit 1; print "Known one-cell heat invariant passed"}'
for mode in blocking overlap; do
  mpirun --bind-to core -np 6 cs1_heat_mpi/parallel 23 51 50 2d "$mode" verify
done
for mode in atomic private; do
  cs2_kmer_openmp/parallel 100003 12 42 uniform "$mode" verify
done
cs3_montecarlo/parallel 1000003 42 good verify
cs3_montecarlo/serial 1000000 42 | awk '
  {for(i=1;i<=NF;i++){split($i,a,"=");v[a[1]]=a[2]}}
  END {x=v["black_scholes"]-10.450583572185565; if(x>1e-10 || x< -1e-10) exit 1; print "Black-Scholes known-value invariant passed"}'
for mode in pull push; do
  cs4_pagerank_openmp/parallel 65536 1000000 42 "$mode" dynamic 64 verify
done
cs4_pagerank_openmp/serial 128 0 42 | awk '
  {for(i=1;i<=NF;i++){split($i,a,"=");v[a[1]]=a[2]}}
  END {if(v["iterations"]!=1 || v["residual"]!=0 || v["rank_sum"]!=1) exit 1; print "All-dangling uniform fixed-point invariant passed"}'
for distribution in uniform skew; do
  for n in 1 100003; do
    mpirun --bind-to core -np 6 cs5_samplesort_mpi/parallel "$n" "$distribution" 8 42 verify
  done
done
echo 'Final six-worker comparisons and known-value invariants passed.'
