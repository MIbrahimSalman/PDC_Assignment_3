include ../common.mk
all: serial parallel pingpong
serial: serial.cpp sort.hpp ../common.hpp
	$(CXX) $(CXXFLAGS) serial.cpp -o $@
parallel: parallel.cpp sort.hpp ../common.hpp
	$(MPICXX) $(CXXFLAGS) parallel.cpp -o $@
pingpong: pingpong.cpp ../common.hpp
	$(MPICXX) $(CXXFLAGS) pingpong.cpp -o $@
.DEFAULT_GOAL := all
