#pragma once
#include "../common.hpp"
#include <fstream>
#include <map>
#ifdef PARALLEL_BUILD
#include <omp.h>
#endif
struct CSR { std::vector<std::size_t> offsets; std::vector<std::uint32_t> neighbors; };
struct Graph { std::uint32_t n; std::size_t edges; CSR outgoing, incoming; };
inline Graph generate_graph(std::uint32_t n, std::size_t edges, std::uint64_t seed) {
    Graph g{n,edges,{}, {}};
    std::vector<std::pair<std::uint32_t,std::uint32_t>> pairs;
    pairs.reserve(edges);
    g.outgoing.offsets.resize(static_cast<std::size_t>(n)+1);
    g.incoming.offsets.resize(static_cast<std::size_t>(n)+1);
    auto state = seed;
    for (std::size_t e = 0; e < edges; ++e) {
        std::uint32_t u = 0, v = 0;
        for (std::uint32_t bit = n/2; bit; bit >>= 1) {
            auto draw = mix64(state++) % 100;
            if (draw >= 76) u |= bit;
            if ((draw >= 57 && draw < 76) || draw >= 95) v |= bit;
        }
        pairs.emplace_back(u,v);
        ++g.outgoing.offsets[u+1];
        ++g.incoming.offsets[v+1];
    }
    std::partial_sum(g.outgoing.offsets.begin(),g.outgoing.offsets.end(),g.outgoing.offsets.begin());
    std::partial_sum(g.incoming.offsets.begin(),g.incoming.offsets.end(),g.incoming.offsets.begin());
    g.outgoing.neighbors.resize(edges);
    g.incoming.neighbors.resize(edges);
    auto out = g.outgoing.offsets, in = g.incoming.offsets;
    for (auto [u,v] : pairs) { g.outgoing.neighbors[out[u]++] = v; g.incoming.neighbors[in[v]++] = u; }
    return g;
}
struct RankResult {
    std::vector<double> ranks;
    double seconds = 0, imbalance = 1, residual = 0;
    int iterations = 0, workers = 1;
};
inline RankResult pagerank(const Graph& g, bool parallel, bool push) {
    RankResult result;
    result.ranks.assign(g.n, 1.0/g.n);
    std::vector<double> next(g.n), contribution(g.n);
    int capacity = 1;
#ifdef PARALLEL_BUILD
    if (parallel) capacity = omp_get_max_threads();
#else
    (void)parallel;
    (void)push;
#endif
    std::vector<double> work_times(capacity,0);
    double start = wall_time();
    for (int iteration = 0; iteration < 500; ++iteration) {
        double dangling = 0;
#ifdef PARALLEL_BUILD
        #pragma omp parallel for if(parallel) reduction(+:dangling) schedule(static)
#endif
        for (std::uint32_t u = 0; u < g.n; ++u) {
            auto degree = g.outgoing.offsets[u+1]-g.outgoing.offsets[u];
            contribution[u] = degree ? result.ranks[u]/degree : 0;
            if (!degree) dangling += result.ranks[u];
        }
        const double base = (1-0.85)/g.n + 0.85*dangling/g.n;
#ifdef PARALLEL_BUILD
        if (parallel) {
            if (push) {
                #pragma omp parallel for schedule(static)
                for (std::uint32_t v = 0; v < g.n; ++v) next[v] = base;
            }
            #pragma omp parallel
            {
                int id = omp_get_thread_num();
                #pragma omp single
                result.workers = omp_get_num_threads();
                double tick = omp_get_wtime();
                #pragma omp for schedule(runtime) nowait
                for (std::uint32_t v = 0; v < g.n; ++v) {
                    if (push) {
                        const double amount = 0.85*contribution[v];
                        for (auto e = g.outgoing.offsets[v]; e < g.outgoing.offsets[v+1]; ++e) {
                            auto destination = g.outgoing.neighbors[e];
                            #pragma omp atomic update
                            next[destination] += amount;
                        }
                    } else {
                        double sum = 0;
                        for (auto e = g.incoming.offsets[v]; e < g.incoming.offsets[v+1]; ++e)
                            sum += contribution[g.incoming.neighbors[e]];
                        next[v] = base+0.85*sum;
                    }
                }
                // Stop before the team's implicit barrier, so this measures work
                // plus scheduling overhead, rather than time waiting for peers.
                work_times[id] += omp_get_wtime()-tick;
            }
        } else
#endif
        {
            for (std::uint32_t v = 0; v < g.n; ++v) {
                double sum = 0;
                for (auto e = g.incoming.offsets[v]; e < g.incoming.offsets[v+1]; ++e)
                    sum += contribution[g.incoming.neighbors[e]];
                next[v] = base+0.85*sum;
            }
        }
        double residual = 0;
#ifdef PARALLEL_BUILD
        #pragma omp parallel for if(parallel) reduction(+:residual) schedule(static)
#endif
        for (std::uint32_t v = 0; v < g.n; ++v) residual += std::abs(next[v]-result.ranks[v]);
        result.ranks.swap(next);
        result.iterations = iteration+1;
        result.residual = residual;
        if (residual < 1e-10) break;
    }
    result.seconds = wall_time()-start;
    double total = std::accumulate(work_times.begin(),work_times.end(),0.0);
    if (total > 0) result.imbalance = *std::max_element(work_times.begin(),work_times.end())/(total/result.workers);
    return result;
}
inline int run(int argc, char** argv, bool parallel) {
    try {
        if (argc > (parallel ? 9 : 5)) throw std::invalid_argument("usage: serial [vertices] [edges] [seed] [degree.csv]; parallel [vertices] [edges] [seed] [pull|push] [static|dynamic|guided] [chunk] [verify|noverify] [degree.csv]");
        auto vertices = argc > 1 ? integer(argv[1]) : 65536;
        auto edges = argc > 2 ? integer(argv[2]) : 1000000;
        auto seed = argc > 3 ? integer(argv[3]) : 42;
        if (vertices < 2 || vertices > (1u<<24) || (vertices & (vertices-1)) || edges > 1000000000)
            throw std::invalid_argument("vertices must be a power of two in [2,2^24]; edges <= 10^9");
        std::string mode = parallel && argc > 4 ? argv[4] : "pull";
        std::string schedule = parallel && argc > 5 ? argv[5] : "static";
        auto chunk = parallel && argc > 6 ? integer(argv[6]) : 0;
        std::string verify = parallel && argc > 7 ? argv[7] : "noverify";
        std::string degree_file = parallel ? (argc > 8 ? argv[8] : "") : (argc > 4 ? argv[4] : "");
        verification_flag(verify);
        if ((mode != "pull" && mode != "push") || (schedule != "static" && schedule != "dynamic" && schedule != "guided") || chunk > 1000000000)
            throw std::invalid_argument("invalid mode, schedule or chunk");
        if (schedule != "static" && chunk == 0) throw std::invalid_argument("dynamic/guided require a positive chunk");
#ifdef PARALLEL_BUILD
        omp_set_dynamic(0);
        omp_set_schedule(schedule == "static" ? omp_sched_static : (schedule == "dynamic" ? omp_sched_dynamic : omp_sched_guided), static_cast<int>(chunk));
#endif
        double start = wall_time();
        auto graph = generate_graph(static_cast<std::uint32_t>(vertices),edges,seed);
        double generation = wall_time()-start;
        if (!degree_file.empty()) {
            std::map<std::size_t,std::size_t> counts;
            for (std::uint32_t v = 0; v < graph.n; ++v) ++counts[graph.incoming.offsets[v+1]-graph.incoming.offsets[v]];
            std::ofstream output(degree_file);
            if (!output) throw std::runtime_error("cannot write degree CSV");
            output << "degree,vertices\n";
            for (auto [degree,count] : counts) output << degree << ',' << count << '\n';
        }
        auto result = pagerank(graph,parallel,mode == "push");
        double error = 0, l1 = 0;
        if (parallel && verify == "verify") {
            auto reference = pagerank(graph,false,false);
            for (std::size_t v = 0; v < vertices; ++v) {
                double difference = std::abs(reference.ranks[v]-result.ranks[v]);
                error = std::max(error,difference); l1 += difference;
            }
        }
        double mass = std::accumulate(result.ranks.begin(),result.ranks.end(),0.0);
        bool valid = result.residual < 1e-10 && error <= 1e-9 && l1 <= 1e-8 && std::abs(mass-1) < 1e-8;
        std::cout << std::setprecision(17) << "seconds=" << result.seconds << " generation_seconds=" << generation
                  << " vertices=" << vertices << " edges=" << edges << " workers=" << result.workers << " mode=" << mode
                  << " schedule=" << schedule << " chunk=" << chunk << " iterations=" << result.iterations
                  << " residual=" << result.residual << " rank_sum=" << mass << " imbalance=" << result.imbalance
                  << " max_error=" << error << " l1_error=" << l1 << " verified="
                  << (verify == "verify" ? (valid ? "yes" : "FAILED") : "not_requested")
                  << " million_edges_per_second=" << static_cast<double>(edges)*result.iterations/result.seconds/1e6 << '\n';
        return valid ? 0 : 2;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
