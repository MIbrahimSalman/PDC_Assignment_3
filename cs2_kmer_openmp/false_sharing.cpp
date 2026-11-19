#include "kmer.hpp"
#include <omp.h>
struct Packed { volatile std::uint64_t value = 0; };
struct alignas(64) Padded { volatile std::uint64_t value = 0; };
template<class Counter> double run(std::uint64_t iterations) {
    std::vector<Counter> counters(omp_get_max_threads());
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        auto& value = counters[omp_get_thread_num()].value;
        // Volatile forces actual accesses instead of collapsing the loop to one store.
        // Each thread owns a distinct object; this experiment has no shared-counter race.
        for (std::uint64_t i = 0; i < iterations; ++i) value = value + 1;
    }
    double elapsed = omp_get_wtime()-start;
    return elapsed;
}
int main(int argc, char** argv) {
    try {
        if (argc > 3) throw std::invalid_argument("usage: false_sharing [iterations=100000000] [packed|padded]");
        auto iterations = argc > 1 ? number(argv[1]) : 100000000;
        std::string mode = argc > 2 ? argv[2] : "packed";
        if (mode != "packed" && mode != "padded") throw std::invalid_argument("invalid counter layout");
        omp_set_dynamic(0);
        auto elapsed = mode == "packed" ? run<Packed>(iterations) : run<Padded>(iterations);
        std::cout << std::setprecision(17) << "seconds=" << elapsed << " layout=" << mode
                  << " iterations=" << iterations << " cache_line_assumption_bytes=64\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
