#include "sort.hpp"
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int rank,p;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&p);
    try {
        if (argc > 6) throw std::invalid_argument("usage: parallel [N] [uniform|skew] [oversampling=1] [seed=42] [verify|noverify]");
        auto n = argc > 1 ? integer(argv[1]) : 100000000;
        std::string distribution = argc > 2 ? argv[2] : "uniform";
        auto oversampling = argc > 3 ? integer(argv[3]) : 1;
        auto seed = argc > 4 ? integer(argv[4]) : 42;
        std::string verify = argc > 5 ? argv[5] : "noverify";
        check_sort_input(n,distribution); verification_flag(verify);
        if (oversampling < 1 || oversampling > 1000 || static_cast<std::uint64_t>(p)*p*oversampling > 10000000)
            throw std::invalid_argument("oversampling/sample metadata too large");
        if (verify == "verify" && n > 1000000) throw std::invalid_argument("exact serial oracle is limited to N<=1000000; production checks use noverify");
        double gen_start = MPI_Wtime();
        auto keys = generate_keys(n,p,rank,seed,distribution == "skew");
        double generation = MPI_Wtime()-gen_start, gen_max;
        MPI_Reduce(&generation,&gen_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        std::uint64_t before[2] = {0,0}, global_before[2];
        for (auto key : keys) { before[0] += key; before[1] += mix64(key); }
        MPI_Allreduce(before,global_before,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        std::vector<int> send_counts(p),recv_counts(p),send_offsets(p),recv_offsets(p);
        std::vector<std::uint32_t> splitters(p-1);
        double phase[4] = {0,0,0,0};
        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime(), tick = start;
        std::sort(keys.begin(),keys.end());
        phase[0] = MPI_Wtime()-tick; tick = MPI_Wtime();
        int sample_count = static_cast<int>(std::min<std::uint64_t>(keys.size(),oversampling*(p-1)));
        std::vector<std::uint32_t> samples(sample_count);
        for (int i = 0; i < sample_count; ++i) samples[i] = keys[(static_cast<std::size_t>(i)+1)*keys.size()/(sample_count+1)];
        std::vector<int> sample_counts(p),sample_offsets(p);
        MPI_Gather(&sample_count,1,MPI_INT,sample_counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
        int sample_total = 0;
        if (rank == 0) for (int i = 0; i < p; ++i) { sample_offsets[i] = sample_total; sample_total += sample_counts[i]; }
        std::vector<std::uint32_t> gathered(sample_total);
        MPI_Gatherv(samples.data(),sample_count,MPI_UINT32_T,gathered.data(),sample_counts.data(),sample_offsets.data(),MPI_UINT32_T,0,MPI_COMM_WORLD);
        if (rank == 0 && p > 1) {
            if (gathered.empty()) throw std::runtime_error("no splitter samples");
            std::sort(gathered.begin(),gathered.end());
            for (int i = 1; i < p; ++i) splitters[i-1] = gathered[static_cast<std::size_t>(i)*gathered.size()/p];
        }
        MPI_Bcast(splitters.data(),p-1,MPI_UINT32_T,0,MPI_COMM_WORLD);
        auto previous = keys.begin();
        for (int i = 0; i < p; ++i) {
            auto end = i == p-1 ? keys.end() : std::upper_bound(previous,keys.end(),splitters[i]);
            send_offsets[i] = static_cast<int>(previous-keys.begin());
            send_counts[i] = static_cast<int>(end-previous);
            previous = end;
        }
        phase[1] = MPI_Wtime()-tick; tick = MPI_Wtime();
        MPI_Alltoall(send_counts.data(),1,MPI_INT,recv_counts.data(),1,MPI_INT,MPI_COMM_WORLD);
        int received = 0;
        for (int i = 0; i < p; ++i) { recv_offsets[i] = received; received += recv_counts[i]; }
        std::vector<std::uint32_t> output(received);
        MPI_Alltoallv(keys.data(),send_counts.data(),send_offsets.data(),MPI_UINT32_T,
                      output.data(),recv_counts.data(),recv_offsets.data(),MPI_UINT32_T,MPI_COMM_WORLD);
        phase[2] = MPI_Wtime()-tick; tick = MPI_Wtime();
        std::sort(output.begin(),output.end());
        phase[3] = MPI_Wtime()-tick;
        MPI_Barrier(MPI_COMM_WORLD);
        double elapsed = MPI_Wtime()-start, max_elapsed, maximum[4];
        MPI_Reduce(&elapsed,&max_elapsed,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        MPI_Reduce(phase,maximum,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        struct { double time; int rank; } local_pair{elapsed,rank}, slowest;
        MPI_Allreduce(&local_pair,&slowest,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
        double critical[4];
        std::copy(phase,phase+4,critical);
        MPI_Bcast(critical,4,MPI_DOUBLE,slowest.rank,MPI_COMM_WORLD);
        std::uint64_t after[3] = {static_cast<std::uint64_t>(output.size()),0,0}, global_after[3];
        for (auto key : output) { after[1] += key; after[2] += mix64(key); }
        MPI_Allreduce(after,global_after,3,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        int ordered = std::is_sorted(output.begin(),output.end()) ? 1 : 0;
        std::vector<int> sizes(p);
        std::vector<std::uint32_t> minima(p),maxima(p);
        std::uint32_t low = output.empty() ? 0 : output.front(), high = output.empty() ? 0 : output.back();
        MPI_Allgather(&received,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
        MPI_Allgather(&low,1,MPI_UINT32_T,minima.data(),1,MPI_UINT32_T,MPI_COMM_WORLD);
        MPI_Allgather(&high,1,MPI_UINT32_T,maxima.data(),1,MPI_UINT32_T,MPI_COMM_WORLD);
        int last = -1;
        for (int i = 0; i < p; ++i) if (sizes[i]) {
            if (last >= 0 && maxima[last] > minima[i]) ordered = 0;
            last = i;
        }
        if (global_after[0] != n || global_after[1] != global_before[0] || global_after[2] != global_before[1]) ordered = 0;
        if (verify == "verify") {
            std::vector<int> offsets(p);
            for (int i = 1; i < p; ++i) offsets[i] = offsets[i-1]+sizes[i-1];
            std::vector<std::uint32_t> actual(rank == 0 ? n : 0);
            MPI_Gatherv(output.data(),received,MPI_UINT32_T,actual.data(),sizes.data(),offsets.data(),MPI_UINT32_T,0,MPI_COMM_WORLD);
            if (rank == 0) {
                auto expected = all_keys(n,p,seed,distribution == "skew");
                std::sort(expected.begin(),expected.end());
                if (actual != expected) ordered = 0;
            }
        }
        int valid;
        MPI_Allreduce(&ordered,&valid,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << std::setprecision(17) << "seconds=" << max_elapsed << " generation_seconds=" << gen_max
                      << " local_sort_seconds=" << critical[0] << " sampling_seconds=" << critical[1]
                      << " exchange_seconds=" << critical[2] << " final_sort_seconds=" << critical[3]
                      << " max_exchange_seconds=" << maximum[2] << " N=" << n << " workers=" << p
                      << " distribution=" << distribution << " oversampling=" << oversampling
                      << " imbalance=" << static_cast<double>(*std::max_element(sizes.begin(),sizes.end()))*p/n
                      << " verified=" << (valid ? "yes" : "FAILED") << " exact_oracle=" << (verify == "verify" ? "yes" : "no") << " bucket_sizes=";
            for (int i = 0; i < p; ++i) std::cout << (i ? "," : "") << sizes[i];
            std::cout << '\n';
        }
        MPI_Finalize(); return valid ? 0 : 2;
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD,1); return 1;
    }
}
