#pragma once
#include "../common.hpp"
#ifdef PARALLEL_BUILD
#include <omp.h>
#endif
constexpr std::uint64_t block_size = 65536;
struct Moments { double sum = 0, squares = 0; };
struct Estimate { double price, standard_error, seconds; int workers; };
inline double closed_form() {
    const double d1 = (0.05 + 0.2*0.2/2)/0.2;
    const double d2 = d1-0.2;
    auto cdf = [](double x) { return 0.5*std::erfc(-x/std::sqrt(2.0)); };
    return 100*cdf(d1)-100*std::exp(-0.05)*cdf(d2);
}
inline Moments sample_block(std::uint64_t start, std::uint64_t end, std::uint64_t seed,
                            std::uint64_t repeat_period) {
    Moments result;
    const auto key = mix64(seed);
    for (auto i = start; i < end; ++i) {
        // Disjoint counter ranges act as a skip-ahead stream assignment. With the
        // bad-stream demonstration enabled, workers intentionally reuse counters.
        auto index = repeat_period ? i % repeat_period : i;
        double u = uniform_open(mix64(key + 2*index));
        double v = uniform_open(mix64(key + 2*index+1));
        double z = std::sqrt(-2*std::log(u))*std::cos(6.2831853071795864769*v);
        double payoff = std::exp(-0.05)*std::max(100*std::exp(0.03+0.2*z)-100, 0.0);
        result.sum += payoff;
        result.squares += payoff*payoff;
    }
    return result;
}
inline Estimate estimate(std::uint64_t samples, std::uint64_t seed, bool parallel, bool bad) {
    std::uint64_t blocks = (samples+block_size-1)/block_size;
    std::vector<Moments> partial(blocks);
    int workers = 1;
    double start = wall_time();
#ifdef PARALLEL_BUILD
    if (parallel) {
        #pragma omp parallel
        {
            #pragma omp single
            workers = omp_get_num_threads();
            // Fixed blocks + deterministic final reduction preserve summation order
            // even when the team size changes. No shared mutable RNG exists.
            #pragma omp for schedule(static)
            for (std::uint64_t b = 0; b < blocks; ++b)
                partial[b] = sample_block(b*block_size, std::min(samples,(b+1)*block_size), seed,
                                         bad ? (samples+workers-1)/workers : 0);
        }
    } else
#else
    (void)parallel;
    (void)bad;
#endif
    {
        for (std::uint64_t b = 0; b < blocks; ++b)
            partial[b] = sample_block(b*block_size, std::min(samples,(b+1)*block_size), seed, 0);
    }
    // Ordered reduction of block partials is intentional: a floating-point OpenMP
    // reduction may choose a different tree between otherwise identical runs.
    long double sum = 0, squares = 0;
    for (const auto& block : partial) { sum += block.sum; squares += block.squares; }
    double mean = static_cast<double>(sum/samples);
    double variance = static_cast<double>((squares-sum*sum/samples)/(samples-1));
    return {mean, std::sqrt(std::max(variance,0.0)/samples), wall_time()-start, workers};
}
inline int run(int argc, char** argv, bool parallel) {
    try {
        if (argc > (parallel ? 5 : 3)) throw std::invalid_argument("usage: program [M=100000000] [seed=42] [good|same_stream] [verify|noverify]");
        auto samples = argc > 1 ? integer(argv[1]) : 100000000;
        auto seed = argc > 2 ? integer(argv[2]) : 42;
        std::string mode = argc > 3 ? argv[3] : "good";
        std::string verify = argc > 4 ? argv[4] : "noverify";
        verification_flag(verify);
        if (samples < 2 || samples > 1000000000000ULL) throw std::invalid_argument("M must be between 2 and 10^12");
        if (mode != "good" && mode != "same_stream") throw std::invalid_argument("invalid stream mode");
#ifdef PARALLEL_BUILD
        omp_set_dynamic(0);
#endif
        auto result = estimate(samples, seed, parallel, mode == "same_stream");
        double error = 0, se_error = 0;
        if (parallel && verify == "verify") {
            auto reference = estimate(samples,seed,false,false);
            error = std::abs(result.price-reference.price);
            se_error = std::abs(result.standard_error-reference.standard_error);
        }
        bool valid = std::isfinite(result.price) && std::isfinite(result.standard_error)
                     && error <= 1e-12 && se_error <= 1e-12;
        std::cout << std::setprecision(17) << "seconds=" << result.seconds << " price=" << result.price
                  << " standard_error=" << result.standard_error << " black_scholes=" << closed_form()
                  << " z_score=" << (result.price-closed_form())/result.standard_error
                  << " samples=" << samples << " workers=" << result.workers << " mode=" << mode
                  << " max_error=" << error << " se_error=" << se_error << " verified="
                  << (verify == "verify" ? (valid ? "yes" : "FAILED") : "not_requested") << '\n';
        return valid ? 0 : 2;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
