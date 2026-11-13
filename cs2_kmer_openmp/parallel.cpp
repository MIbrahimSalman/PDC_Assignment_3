#include "kmer.hpp"
#include <omp.h>
int main(int argc, char** argv) {
    try {
        if (argc > 7) throw std::invalid_argument("usage: parallel [L] [k] [seed] [uniform|skew] [race|atomic|private] [verify|noverify]");
        auto x = arguments(argc, argv);
        std::string mode = argc > 5 ? argv[5] : "private";
        std::string verify = argc > 6 ? argv[6] : "verify";
        if ((mode != "race" && mode != "atomic" && mode != "private") ||
            (verify != "verify" && verify != "noverify")) throw std::invalid_argument("invalid mode or verification flag");
        omp_set_dynamic(0);
        double start = seconds();
        auto s = generate(x);
        double generation = seconds()-start;
        Histogram h(x.bins(), 0);
        std::vector<Histogram> private_hist;
        int workers = 1;
        start = seconds();
        if (mode == "private") {
            // Include private allocation, zeroing, and merge in the parallel time.
            private_hist.resize(omp_get_max_threads());
            for (auto& histogram : private_hist) histogram.resize(x.bins(), 0);
        }
        if (mode == "race") {
            // Deliberately incorrect assignment demonstration. Data races are undefined
            // behavior; even an apparently correct total does NOT establish safety.
            const auto blocks = static_cast<std::uint64_t>(omp_get_max_threads());
            #pragma omp parallel
            {
                #pragma omp single
                workers = omp_get_num_threads();
                #pragma omp for schedule(static)
                for (std::uint64_t block = 0; block < blocks; ++block) {
                    const auto q = x.windows()/blocks, r = x.windows()%blocks;
                    auto begin = block*q + std::min(block, r);
                    auto end = begin+q+(block < r);
                    count_range(s, x.k, begin, end, [&](std::uint32_t code) { ++h[code]; });
                }
            }
        } else {
            #pragma omp parallel
            {
                int id = omp_get_thread_num(), p = omp_get_num_threads();
                #pragma omp single
                workers = p;
                const auto q = x.windows()/p, r = x.windows()%p;
                auto begin = id*q+std::min(static_cast<std::uint64_t>(id), r);
                auto end = begin+q+(static_cast<std::uint64_t>(id) < r);
                if (mode == "atomic") {
                    count_range(s, x.k, begin, end, [&](std::uint32_t code) {
                        #pragma omp atomic update
                        ++h[code];
                    });
                } else {
                    auto& local = private_hist[id];
                    count_range(s, x.k, begin, end, [&](std::uint32_t code) { ++local[code]; });
                }
            }
            if (mode == "private") {
                #pragma omp parallel for schedule(static)
                for (std::size_t bin = 0; bin < x.bins(); ++bin)
                    for (int id = 0; id < workers; ++id) h[bin] += private_hist[id][bin];
            }
        }
        double elapsed = seconds()-start;
        auto total = std::accumulate(h.begin(), h.end(), std::uint64_t{0});
        bool correct = total == x.windows();
        std::size_t mismatches = 0;
        if (verify == "verify") {
            Histogram expected(x.bins(), 0);
            serial_count(s, x, expected);
            for (std::size_t i = 0; i < h.size(); ++i) mismatches += h[i] != expected[i];
            correct = correct && mismatches == 0;
        }
        std::cout << std::setprecision(17) << "seconds=" << elapsed << " generation_seconds=" << generation
                  << " total=" << total << " expected=" << x.windows() << " mismatched_bins=" << mismatches
                  << " verified=" << (verify == "verify" ? (correct ? "yes" : "FAILED") : "not_requested")
                  << " workers=" << workers << " mode=" << mode << " L=" << x.length << " k=" << x.k
                  << " private_bytes=" << (private_hist.size()*x.bins()*sizeof(std::uint64_t))
                  << " million_kmers_per_second=" << x.windows()/elapsed/1e6 << '\n';
        // Version A's wrong answer is observable; do not disguise it as a passing check.
        return correct ? 0 : 2;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
