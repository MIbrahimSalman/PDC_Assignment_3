#include "heat.hpp"
int main(int argc, char** argv) {
    try {
        if (argc > 4) throw std::invalid_argument("usage: serial [N=1024] [steps=500] [K=50]");
        int n = argc > 1 ? positive(argv[1]) : 1024;
        int t = argc > 2 ? positive(argv[2]) : 500;
        int k = argc > 3 ? positive(argv[3]) : 50;
        if (n < 3) throw std::invalid_argument("N must be at least 3");
        auto a = grid(n), b = a;
        double start = seconds();
        double delta = solve_serial(a, b, n, t, k);
        double elapsed = seconds() - start;
        double sum = 0.0;
        for (double v : a) sum += v;
        std::cout << std::setprecision(17) << "seconds=" << elapsed << " max_change=" << delta
                  << " checksum=" << sum << " N=" << n << " steps=" << t << " K=" << k << '\n';
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
