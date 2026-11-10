#include <iostream>
#include <omp.h>
int main() { std::cout << "omp_get_max_threads()=" << omp_get_max_threads() << '\n'; }
