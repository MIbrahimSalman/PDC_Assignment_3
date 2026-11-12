include ../common.mk
all: serial parallel
serial: serial.cpp montecarlo.hpp ../common.hpp
	$(CXX) $(CXXFLAGS) serial.cpp -o $@
parallel: parallel.cpp montecarlo.hpp ../common.hpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) parallel.cpp -o $@
.DEFAULT_GOAL := all
