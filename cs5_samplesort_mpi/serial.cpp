#include "sort.hpp"
int main(int argc, char** argv) {
    try {
        if (argc > 5) throw std::invalid_argument("usage: serial [N=100000000] [shards=1] [uniform|skew] [seed=42]");
        auto n = argc > 1 ? integer(argv[1]) : 100000000;
        auto shards = argc > 2 ? integer(argv[2]) : 1;
        std::string distribution = argc > 3 ? argv[3] : "uniform";
        auto seed = argc > 4 ? integer(argv[4]) : 42;
        check_sort_input(n,distribution);
        if (shards < 1 || shards > 1000000) throw std::invalid_argument("invalid shard count");
        double start = wall_time();
        auto keys = all_keys(n,static_cast<int>(shards),seed,distribution == "skew");
        double generation = wall_time()-start;
        start = wall_time();
        std::sort(keys.begin(),keys.end());
        double elapsed = wall_time()-start;
        bool correct = std::is_sorted(keys.begin(),keys.end()) && keys.size() == n;
        std::cout << std::setprecision(17) << "seconds=" << elapsed << " generation_seconds=" << generation
                  << " N=" << n << " shards=" << shards << " distribution=" << distribution
                  << " verified=" << (correct ? "yes" : "FAILED") << '\n';
        return correct ? 0 : 2;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
