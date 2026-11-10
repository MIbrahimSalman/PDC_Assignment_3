CXX = g++
MPICXX = mpic++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Wpedantic
OMPFLAGS = -fopenmp

.PHONY: all clean
clean:
	rm -f serial parallel false_sharing hello pingpong *.exe *.o
