#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
g++ -O3 -std=c++17 -fopenmp tools/hello.cpp -o tools/hello
OMP_NUM_THREADS=2 tools/hello
mpicc --version | head -n 1
mpirun --version | head -n 2
# Small, uneven decompositions and one-cell-wide local blocks.
for p in 1 2 4; do
  for n in 6 19; do
    for decomposition in 1d 2d; do
      for mode in blocking overlap; do
        for k in 1 50; do
          mpirun -np "$p" cs1_heat_mpi/parallel "$n" 17 "$k" "$decomposition" "$mode" verify
        done
      done
    done
  done
done
mpirun -np 1 cs1_heat_mpi/parallel 3 1 1 2d overlap verify
for p in 1 2 4; do
  for specification in '8 8' '11 8' '10003 1' '100003 8' '100003 12'; do
    read -r length k <<< "$specification"
    for distribution in uniform skew; do
      for mode in atomic private; do
        OMP_NUM_THREADS="$p" cs2_kmer_openmp/parallel "$length" "$k" 42 "$distribution" "$mode" verify
      done
    done
  done
done
echo 'All correct variants passed. The intentional race is excluded from pass/fail tests.'
