#include "block_matmul.h"

struct Config {
	/* MPI Files */
	MPI_File A_file, B_file, C_file;
	char *outfile;

	/* MPI Datatypes for matrix blocks */
	MPI_Datatype block;

	/* Matrix data */
	double *A, *A_tmp, *B, *C;

	/* Cart communicators */
	MPI_Comm grid_comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;

	/* Cart communicator dim and ranks */
	int dim[2], coords[2];
	int world_rank, world_size, grid_rank;
	int row_rank, row_size, col_rank, col_size;

	/* Full matrix dim */
	int A_dims[2];
	int B_dims[2];
	int C_dims[2];
	int matrix_size;

	/* Process local matrix dim */
	int local_dims[2];
	int local_size;
};

struct Config config;

void init_matmul(char *A_file, char *B_file, char *outfile)
{
	/* Copy output file name to configuration */
	//config.outfile = outfile;
	//int world_rank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	/* Get matrix size header */
	printf("%d", 1);
	/*if(world_rank == 0){
		MPI_File fh;
		MPI_Offset offset;

		MPI_File_open(MPI_COMM_SELF, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
		MPI_File_read_at(fh, 0, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);
		printf("%d", config.A_dims[0]);
	}*/
	

	/* Broadcast global matrix sizes */

	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */

	/* Verify dim of A and B matches for matul and both are square*/

	/* Create Cart communicator for NxN processes */

	/* Sub div cart communicator to N row communicator */

	/* Sub div cart communicator to N col communicator */

	/* Setup sizes of full matrices */

	/* Setup sizes of local matrix tiles */

	/* Create subarray datatype for local matrix tile */

	/* Create data array to load actual block matrix data */

	/* Set fileview of process to respective matrix block */

	/* Collective read blocks from files */

	/* Close data source files */
}

void cleanup_matmul()
{
	/* Rank zero writes header specifying dim of result matrix C */

	/* Set fileview of process to respective matrix block with header offset */

	/* Collective write and close file */

	/* Cleanup */
}

void compute_fox()
{

	/* Compute source and target for verticle shift of B blocks */
	int i;
	for (i = 0; i < config.dim[0]; i++) {
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */

		/* dgemm with blocks */
		
		/* Shfting block B upwards and receive from process below */

	}
}
