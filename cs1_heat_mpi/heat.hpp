#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// N includes fixed boundary cells. alpha <= 1/4 makes this explicit scheme stable.
constexpr double alpha = 0.2;
inline double seconds() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
inline int positive(const char* text) {
    std::size_t used = 0;
    long long value = std::stoll(text, &used);
    if (used != std::string(text).size() || value < 1 || value > 1000000000)
        throw std::invalid_argument("expected a positive integer <= 1000000000");
    return static_cast<int>(value);
}
inline double initial(int row, int /*column*/) { return row == 0 ? 100.0 : 0.0; }
inline std::vector<double> grid(int n) {
    std::vector<double> result(static_cast<std::size_t>(n) * n, 0.0);
    std::fill_n(result.begin(), n, 100.0);
    return result;
}
inline double update(double center, double north, double south, double west, double east) {
    return center + alpha * (north + south + west + east - 4.0 * center);
}
inline double solve_serial(std::vector<double>& a, std::vector<double>& b, int n, int steps, int k) {
    double change = 0.0;
    for (int t = 1; t <= steps; ++t) {
        const bool check = t % k == 0 || t == steps;
        double local = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                const auto x = static_cast<std::size_t>(i) * n + j;
                b[x] = update(a[x], a[x-n], a[x+n], a[x-1], a[x+1]);
                if (check) local = std::max(local, std::abs(b[x] - a[x]));
            }
        }
        a.swap(b);
        if (check) change = local;
    }
    return change;
}
inline int extent(int total, int parts, int coord) { return total / parts + (coord < total % parts); }
inline int offset(int total, int parts, int coord) { return coord * (total / parts) + std::min(coord, total % parts); }
