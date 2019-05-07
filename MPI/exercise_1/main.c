
#include <mpi.h>



int main(int argc, char **argv)
{
    int rank, num_ranks, prev_rank, next_rank, recv_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    prev_rank = (rank + (num_ranks -1)) % num_ranks;
    next_rank = (rank + 1) % num_ranks;

    MPI_Send(&rank, 1, MPI_INT, next_rank, rank, MPI_COMM_WORLD);
    MPI_Recv(&recv_rank, 1, MPI_INT, prev_rank, prev_rank, MPI_COMM_WORLD, ...);

    printf("[%d] I have received %d from rank=%d\n", rank, recv_rank, prev_rank);
    MPI_Finalize();
		
}


