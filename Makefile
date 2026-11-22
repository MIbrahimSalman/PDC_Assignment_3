CASES = cs1_heat_mpi cs2_kmer_openmp cs3_montecarlo cs4_pagerank_openmp cs5_samplesort_mpi
.PHONY: all clean check
all:
	@for dir in $(CASES); do $(MAKE) -C $$dir || exit $$?; done
check: all
	bash tools/check.sh
	bash tools/check_remaining.sh
clean:
	@for dir in $(CASES); do $(MAKE) -C $$dir clean || exit $$?; done
	rm -f tools/hello
