#pragma once
#include "../common.hpp"
inline std::vector<std::uint32_t> generate_keys(std::uint64_t n, int p, int rank, std::uint64_t seed, bool skew) {
    auto count = partition_begin(n,p,rank+1)-partition_begin(n,p,rank);
    std::vector<std::uint32_t> keys(count);
    auto state = mix64(seed+rank);
    for (auto& key : keys) {
        auto bits = mix64(state++);
        key = skew && bits % 100 < 90 ? 0 : static_cast<std::uint32_t>(bits >> 32);
    }
    return keys;
}
inline void check_sort_input(std::uint64_t n, const std::string& distribution) {
    if (n < 1 || n > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("N must be in [1,INT_MAX] for MPI_Alltoallv int displacements");
    if (distribution != "uniform" && distribution != "skew") throw std::invalid_argument("expected uniform or skew");
}
inline std::vector<std::uint32_t> all_keys(std::uint64_t n, int p, std::uint64_t seed, bool skew) {
    std::vector<std::uint32_t> keys;
    keys.reserve(n);
    for (int rank = 0; rank < p; ++rank) {
        auto shard = generate_keys(n,p,rank,seed,skew);
        keys.insert(keys.end(),shard.begin(),shard.end());
    }
    return keys;
}
