#!/usr/bin/env bash
# Real measurements only. Run from an otherwise idle machine. No concurrent jobs.
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
mkdir -p results
make > results/build.log 2>&1
{
  date -Is
  uname -a
  g++ --version
  mpicc --version
  mpirun --version
  lscpu
  free -b
  printenv OMP_PROC_BIND OMP_PLACES OMP_DYNAMIC
} > results/environment.txt
# Refuse to mix measurements from different source versions on resume.
find . -type f \( -name '*.cpp' -o -name '*.hpp' -o -name '*.mk' -o -name Makefile \) -print0 | sort -z | xargs -0 sha256sum > results/current_sources.sha256
if [[ -f results/sources.sha256 ]]; then
  cmp results/sources.sha256 results/current_sources.sha256 || { echo 'Source changed: archive results before restarting.'; exit 1; }
else
  cp results/current_sources.sha256 results/sources.sha256
fi
touch results/raw.tsv
run() {
  local label="$1" p="$2" repetition="$3"
  shift 3
  local key
  printf -v key '%s\t%s\t%s\t' "$label" "$p" "$repetition"
  if grep -Fq "$key" results/raw.tsv; then return; fi
  echo "Running $label p=$p repeat=$repetition"
  local output code=0
  output=$(OMP_NUM_THREADS="$p" "$@" 2>>results/stderr.log) || code=$?
  if [[ "$label" == race ]]; then
    [[ "$code" == 0 || "$code" == 2 ]] || return "$code"
  elif [[ "$code" != 0 ]]; then
    echo "$output" >&2
    return "$code"
  fi
  [[ "$output" == *seconds=* ]] || { echo "No measurement: $output"; return 1; }
  printf '%s%s\n' "$key" "$output" >> results/raw.tsv
}
mpi() { mpirun --bind-to core --map-by core -np "$@"; }
# Repetitions are outermost to reduce systematic order/thermal bias.
for r in 1 2 3; do
  for n in 1024 4096; do
    run "heat${n}_serial" 1 "$r" cs1_heat_mpi/serial "$n" 500 50
    for p in 1 2 4 6; do
      run "heat${n}_strong" "$p" "$r" mpi "$p" cs1_heat_mpi/parallel "$n" 500 50 2d overlap noverify
    done
  done
  for d in 1d 2d; do
    for mode in blocking overlap; do
      for k in 1 50; do
        run "heat_compare_${d}_${mode}_k${k}" 4 "$r" mpi 4 cs1_heat_mpi/parallel 1024 500 "$k" "$d" "$mode" noverify
      done
    done
  done
  for p in 1 2 4 6; do
    n=$(awk -v p="$p" 'BEGIN {printf "%d", 2+int(512*sqrt(p)+0.5)}')
    run heat_weak "$p" "$r" mpi "$p" cs1_heat_mpi/parallel "$n" 500 50 2d overlap noverify
  done
  for k in 8 12; do
    run "kmer${k}_serial" 1 "$r" cs2_kmer_openmp/serial 100000000 "$k" 42 uniform
    for p in 1 2 4 6; do
      for mode in atomic private; do
        run "kmer${k}_${mode}" "$p" "$r" cs2_kmer_openmp/parallel 100000000 "$k" 42 uniform "$mode" noverify
      done
    done
  done
  run race 4 "$r" cs2_kmer_openmp/parallel 1000000 8 42 skew race verify
  for p in 1 2 4 6; do
    for mode in packed padded; do
      run "false_sharing_${mode}" "$p" "$r" cs2_kmer_openmp/false_sharing 100000000 "$mode"
    done
  done
  run mc_serial 1 "$r" cs3_montecarlo/serial 100000000 42
  for p in 1 2 4 6 8; do
    run mc_good "$p" "$r" cs3_montecarlo/parallel 100000000 42 good noverify
  done
  run mc_same_stream 8 "$r" cs3_montecarlo/parallel 100000000 42 same_stream noverify
  run pagerank_serial 1 "$r" cs4_pagerank_openmp/serial 65536 1000000 42 results/degrees.csv
  for p in 1 2 4 6; do
    run pagerank_strong "$p" "$r" cs4_pagerank_openmp/parallel 65536 1000000 42 pull static 0 noverify
  done
  for mode in pull push; do
    for schedule in static dynamic guided; do
      for chunk in 1 64 256; do
        run "pagerank_${mode}_${schedule}_${chunk}" 4 "$r" cs4_pagerank_openmp/parallel 65536 1000000 42 "$mode" "$schedule" "$chunk" noverify
      done
    done
  done
  for distribution in uniform skew; do
    for p in 1 2 4 6; do
      # Generate identical shards in serial to ensure a fair data comparison.
      run "sort_${distribution}_serial" "$p" "$r" cs5_samplesort_mpi/serial 100000000 "$p" "$distribution" 42
      run "sort_${distribution}_strong" "$p" "$r" mpi "$p" cs5_samplesort_mpi/parallel 100000000 "$distribution" 1 42 noverify
    done
    run "sort_${distribution}_oversample8" 4 "$r" mpi 4 cs5_samplesort_mpi/parallel 100000000 "$distribution" 8 42 noverify
  done
  mpi 2 cs5_samplesort_mpi/pingpong > "results/pingpong_${r}.csv"
done
# Calibrate equal target runtime from median observed throughput, then measure it.
# This is approximate fixed-time work: actual times and mismatch are reported.
for p in 1 8; do
  t=$(awk -F '\t' -v p="$p" '$1=="mc_good" && $2==p {split($4,a," "); split(a[1],b,"="); print b[2]}' results/raw.tsv | sort -n | sed -n '2p')
  samples=$(awk -v t="$t" 'BEGIN {printf "%.0f", 3*100000000/t}')
  for r in 1 2 3; do
    run mc_fixed_time "$p" "$r" cs3_montecarlo/parallel "$samples" 42 good noverify
  done
done
## Additional final-size fixed-time comparison (exploratory 3-second data retained).
bash tools/fixed_time.sh
echo 'All planned benchmark runs completed.'
