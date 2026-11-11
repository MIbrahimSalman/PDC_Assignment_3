#include "kmer.hpp"
int main(int argc, char** argv) {
    try {
        if (argc > 5) throw std::invalid_argument("usage: serial [L] [k] [seed] [uniform|skew]");
        auto x = arguments(argc, argv);
        double start = seconds();
        auto s = generate(x);
        double generation = seconds()-start;
        Histogram h(x.bins(), 0);
        start = seconds();
        serial_count(s, x, h);
        double elapsed = seconds()-start;
        auto total = std::accumulate(h.begin(), h.end(), std::uint64_t{0});
        std::cout << std::setprecision(17) << "seconds=" << elapsed << " generation_seconds=" << generation
                  << " total=" << total << " expected=" << x.windows() << " L=" << x.length << " k=" << x.k
                  << " million_kmers_per_second=" << x.windows()/elapsed/1e6 << '\n';
        return total == x.windows() ? 0 : 2;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
