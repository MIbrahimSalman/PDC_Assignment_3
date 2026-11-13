#include "heat.hpp"
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
#include <array>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    try {
        if (argc > 7) throw std::invalid_argument("usage: parallel [N] [steps] [K] [1d|2d] [blocking|overlap] [verify|noverify]");
        int n = argc > 1 ? positive(argv[1]) : 1024;
        int steps = argc > 2 ? positive(argv[2]) : 500;
        int k = argc > 3 ? positive(argv[3]) : 50;
        std::string decomposition = argc > 4 ? argv[4] : "2d";
        std::string mode = argc > 5 ? argv[5] : "overlap";
        std::string verification = argc > 6 ? argv[6] : "noverify";
        if (n < 3 || (decomposition != "1d" && decomposition != "2d") ||
            (mode != "blocking" && mode != "overlap") ||
            (verification != "verify" && verification != "noverify"))
            throw std::invalid_argument("invalid grid, decomposition, exchange mode, or verification flag");
        int dims[2] = {0, 0}, periods[2] = {0, 0};
        if (decomposition == "1d") { dims[0] = p; dims[1] = 1; }
        else MPI_Dims_create(p, 2, dims);
        if (dims[0] > n-2 || dims[1] > n-2) throw std::invalid_argument("too many ranks for interior grid");
        MPI_Comm cart;
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);
        int coords[2], north, south, west, east;
        MPI_Cart_coords(cart, rank, 2, coords);
        MPI_Cart_shift(cart, 0, 1, &north, &south);
        MPI_Cart_shift(cart, 1, 1, &west, &east);
        int rows = extent(n-2, dims[0], coords[0]);
        int cols = extent(n-2, dims[1], coords[1]);
        int row0 = 1 + offset(n-2, dims[0], coords[0]);
        int col0 = 1 + offset(n-2, dims[1], coords[1]);
        const std::size_t stride = static_cast<std::size_t>(cols) + 2;
        std::vector<double> a((static_cast<std::size_t>(rows)+2)*stride), b(a.size());
        auto index = [stride](int i, int j) { return static_cast<std::size_t>(i)*stride+j; };
        for (int i = 0; i <= rows+1; ++i)
            for (int j = 0; j <= cols+1; ++j) a[index(i,j)] = initial(row0+i-1, col0+j-1);
        b = a;
        MPI_Datatype column;
        MPI_Type_vector(rows, 1, cols+2, MPI_DOUBLE, &column);
        MPI_Type_commit(&column);
        double compute = 0, exchange = 0, reductions = 0, delta = 0;
        MPI_Barrier(cart);
        double start = MPI_Wtime();
        for (int t = 1; t <= steps; ++t) {
            const bool check = t % k == 0 || t == steps;
            double local = 0;
            auto cell = [&](int i, int j) {
                auto x = index(i,j);
                b[x] = update(a[x], a[x-stride], a[x+stride], a[x-1], a[x+1]);
                if (check) local = std::max(local, std::abs(b[x]-a[x]));
            };
            double tick = MPI_Wtime();
            if (mode == "blocking") {
                MPI_Sendrecv(&a[index(1,1)], cols, MPI_DOUBLE, north, 0,
                             &a[index(rows+1,1)], cols, MPI_DOUBLE, south, 0, cart, MPI_STATUS_IGNORE);
                MPI_Sendrecv(&a[index(rows,1)], cols, MPI_DOUBLE, south, 1,
                             &a[index(0,1)], cols, MPI_DOUBLE, north, 1, cart, MPI_STATUS_IGNORE);
                MPI_Sendrecv(&a[index(1,1)], 1, column, west, 2,
                             &a[index(1,cols+1)], 1, column, east, 2, cart, MPI_STATUS_IGNORE);
                MPI_Sendrecv(&a[index(1,cols)], 1, column, east, 3,
                             &a[index(1,0)], 1, column, west, 3, cart, MPI_STATUS_IGNORE);
                exchange += MPI_Wtime()-tick;
                tick = MPI_Wtime();
                for (int i = 1; i <= rows; ++i) for (int j = 1; j <= cols; ++j) cell(i,j);
                compute += MPI_Wtime()-tick;
            } else {
                std::array<MPI_Request,8> req;
                MPI_Irecv(&a[index(rows+1,1)], cols, MPI_DOUBLE, south, 0, cart, &req[0]);
                MPI_Irecv(&a[index(0,1)], cols, MPI_DOUBLE, north, 1, cart, &req[1]);
                MPI_Irecv(&a[index(1,cols+1)], 1, column, east, 2, cart, &req[2]);
                MPI_Irecv(&a[index(1,0)], 1, column, west, 3, cart, &req[3]);
                MPI_Isend(&a[index(1,1)], cols, MPI_DOUBLE, north, 0, cart, &req[4]);
                MPI_Isend(&a[index(rows,1)], cols, MPI_DOUBLE, south, 1, cart, &req[5]);
                MPI_Isend(&a[index(1,1)], 1, column, west, 2, cart, &req[6]);
                MPI_Isend(&a[index(1,cols)], 1, column, east, 3, cart, &req[7]);
                exchange += MPI_Wtime()-tick;
                tick = MPI_Wtime();
                for (int i = 2; i < rows; ++i) for (int j = 2; j < cols; ++j) cell(i,j);
                compute += MPI_Wtime()-tick;
                tick = MPI_Wtime();
                MPI_Waitall(8, req.data(), MPI_STATUSES_IGNORE);
                exchange += MPI_Wtime()-tick;
                tick = MPI_Wtime();
                // Border cells are visited exactly once, even for one-cell-wide blocks.
                for (int i = 1; i <= rows; ++i) {
                    cell(i,1);
                    if (cols > 1) cell(i,cols);
                }
                for (int j = 2; j < cols; ++j) {
                    cell(1,j);
                    if (rows > 1) cell(rows,j);
                }
                compute += MPI_Wtime()-tick;
            }
            if (check) {
                tick = MPI_Wtime();
                MPI_Allreduce(&local, &delta, 1, MPI_DOUBLE, MPI_MAX, cart);
                reductions += MPI_Wtime()-tick;
            }
            a.swap(b);
        }
        MPI_Barrier(cart);
        double timing[4] = {MPI_Wtime()-start, compute, exchange, reductions}, maximum[4];
        MPI_Reduce(timing, maximum, 4, MPI_DOUBLE, MPI_MAX, 0, cart);
        double error = 0, global_error = 0;
        if (verification == "verify") {
            // Debug-only: replicate the serial oracle; never included in timing.
            auto reference = grid(n), scratch = reference;
            solve_serial(reference, scratch, n, steps, k);
            for (int i = 1; i <= rows; ++i) for (int j = 1; j <= cols; ++j)
                error = std::max(error, std::abs(a[index(i,j)]-reference[static_cast<std::size_t>(row0+i-1)*n+col0+j-1]));
            MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, cart);
        }
        if (rank == 0) {
            std::cout << std::setprecision(17) << "seconds=" << maximum[0] << " compute_seconds=" << maximum[1]
                      << " exchange_seconds=" << maximum[2] << " reduction_seconds=" << maximum[3]
                      << " max_change=" << delta << " N=" << n << " steps=" << steps << " K=" << k
                      << " workers=" << p << " decomposition=" << decomposition << " mode=" << mode
                      << " verified=" << (verification == "verify" ? (global_error <= 1e-12 ? "yes" : "FAILED") : "not_requested")
                      << " max_error=" << global_error << '\n';
        }
        MPI_Type_free(&column);
        MPI_Comm_free(&cart);
        MPI_Finalize();
        return global_error > 1e-12 ? 2 : 0;
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
