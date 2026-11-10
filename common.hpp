#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
inline double wall_time() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
inline std::uint64_t integer(const char* s) {
    std::string text(s);
    if (text.empty() || text[0] == '-') throw std::invalid_argument("expected nonnegative integer");
    std::size_t used;
    auto n = std::stoull(text, &used);
    if (used != text.size()) throw std::invalid_argument("invalid integer");
    return n;
}
inline std::uint64_t mix64(std::uint64_t z) {
    z += UINT64_C(0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
inline double uniform_open(std::uint64_t bits) {
    // 52 random bits, strictly inside (0,1), including after rounding.
    return (static_cast<double>(bits >> 12) + 0.5) * 0x1.0p-52;
}
inline std::uint64_t partition_begin(std::uint64_t n, int p, int rank) {
    return n/p*rank + std::min(n%p, static_cast<std::uint64_t>(rank));
}
inline void verification_flag(const std::string& flag) {
    if (flag != "verify" && flag != "noverify") throw std::invalid_argument("expected verify or noverify");
}
