include ../common.mk
all: serial parallel
serial: serial.cpp heat.hpp
	$(CXX) $(CXXFLAGS) serial.cpp -o $@
parallel: parallel.cpp heat.hpp
	$(MPICXX) $(CXXFLAGS) parallel.cpp -o $@
.DEFAULT_GOAL := all
