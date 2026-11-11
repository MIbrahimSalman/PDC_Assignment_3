include ../common.mk
all: serial parallel false_sharing
serial: serial.cpp kmer.hpp
	$(CXX) $(CXXFLAGS) serial.cpp -o $@
parallel: parallel.cpp kmer.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) parallel.cpp -o $@
false_sharing: false_sharing.cpp kmer.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) false_sharing.cpp -o $@
.DEFAULT_GOAL := all
