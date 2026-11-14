include ../common.mk
all: serial parallel
serial: serial.cpp pagerank.hpp ../common.hpp
	$(CXX) $(CXXFLAGS) serial.cpp -o $@
parallel: parallel.cpp pagerank.hpp ../common.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) parallel.cpp -o $@
.DEFAULT_GOAL := all
