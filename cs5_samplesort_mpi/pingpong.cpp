#include "../common.hpp"
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int rank,p; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&p);
    if (p != 2) { if (!rank) std::cerr << "pingpong requires exactly 2 ranks\n"; MPI_Finalize(); return 1; }
    if (!rank) std::cout << "bytes,one_way_seconds\n";
    for (int size : {0,8,64,512,4096,32768,262144,1048576,4194304}) {
        std::vector<char> message(std::max(size,1),0);
        int repetitions = size < 262144 ? 1000 : 100;
        for (int round = 0; round < 2; ++round) {
            int count = round == 0 ? 10 : repetitions;
            MPI_Barrier(MPI_COMM_WORLD);
            double start = MPI_Wtime();
            for (int i = 0; i < count; ++i) {
                if (rank == 0) {
                    MPI_Send(message.data(),size,MPI_BYTE,1,0,MPI_COMM_WORLD);
                    MPI_Recv(message.data(),size,MPI_BYTE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(message.data(),size,MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    MPI_Send(message.data(),size,MPI_BYTE,0,0,MPI_COMM_WORLD);
                }
            }
            double elapsed = MPI_Wtime()-start;
            if (rank == 0 && round == 1) std::cout << size << ',' << std::setprecision(17) << elapsed/(2*count) << '\n';
        }
    }
    MPI_Finalize(); return 0;
}
