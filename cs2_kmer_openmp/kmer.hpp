#pragma once
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using Histogram = std::vector<std::uint64_t>;
inline double seconds() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
inline std::uint64_t number(const char* text) {
    std::string s(text);
    if (s.empty() || s[0] == '-') throw std::invalid_argument("expected a nonnegative integer");
    std::size_t used;
    auto x = std::stoull(s, &used);
    if (used != s.size()) throw std::invalid_argument("invalid integer");
    return x;
}
struct Input {
    std::uint64_t length = 100000000, seed = 42;
    int k = 8;
    bool skew = false;
    std::size_t bins() const { return std::size_t{1} << (2*k); }
    std::uint64_t windows() const { return length-static_cast<unsigned>(k)+1; }
};
inline Input arguments(int argc, char** argv) {
    Input x;
    if (argc > 1) x.length = number(argv[1]);
    if (argc > 2) {
        auto k = number(argv[2]);
        if (k < 1 || k > 12) throw std::invalid_argument("k must be between 1 and 12");
        x.k = static_cast<int>(k);
    }
    if (argc > 3) x.seed = number(argv[3]);
    if (argc > 4) {
        std::string distribution(argv[4]);
        if (distribution != "uniform" && distribution != "skew")
            throw std::invalid_argument("distribution must be uniform or skew");
        x.skew = distribution == "skew";
    }
    if (x.length < static_cast<unsigned>(x.k)) throw std::invalid_argument("length must be >= k");
    return x;
}
// SplitMix64: deterministic synthetic input, independent of worker count.
inline std::uint64_t random64(std::uint64_t& state) {
    auto z = (state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
inline std::vector<std::uint8_t> generate(const Input& x) {
    std::vector<std::uint8_t> sequence(x.length);
    auto state = x.seed;
    for (auto& base : sequence) {
        auto r = random64(state);
        // 90% A; the remaining 10% split among C/G/T.
        base = x.skew ? (r % 100 < 90 ? 0 : 1 + (r >> 32) % 3) : r & 3;
    }
    return sequence;
}
template<class Increment>
inline void count_range(const std::vector<std::uint8_t>& s, int k, std::uint64_t begin,
                        std::uint64_t end, Increment increment) {
    if (begin == end) return;
    std::uint32_t code = 0, mask = (std::uint32_t{1} << (2*k))-1;
    for (int j = 0; j < k; ++j) code = (code << 2) | s[begin+j];
    increment(code);
    for (auto i = begin+1; i < end; ++i) {
        code = ((code << 2) | s[i+k-1]) & mask;
        increment(code);
    }
}
inline void serial_count(const std::vector<std::uint8_t>& s, const Input& x, Histogram& h) {
    count_range(s, x.k, 0, x.windows(), [&](std::uint32_t code) { ++h[code]; });
}
