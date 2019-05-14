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

void init_matmuls(char *A_file, char *B_file, char *outfile)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config.world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	config.outfile = outfile;

	int wrap[2];
	wrap[0] = wrap[1] = 1;
	config.dim[0] = 8;
	config.dim[1] = 8;
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, wrap, 1, &config.grid_comm);
	if(config.world_rank == 63){
		int coord[2];
		MPI_Cart_coords(config.grid_comm, config.world_rank, 2, coord);
		//printf("%d\n", coord[0]);
		//printf("%d\n", coord[1]);
		MPI_Comm_rank(config.grid_comm, &config.grid_rank);
		//printf("%d\n", config.grid_rank);
	}
	//MPI_Comm_rank(config.grid_comm, &config.grid_rank);
	config.coords[0] = 0;
	config.coords[1] = 1;
	MPI_Cart_sub(config.grid_comm, config.coords, &config.row_comm);
	MPI_Comm_rank(config.row_comm, &config.row_rank);
	config.coords[0] = 1;
	config.coords[1] = 0;
	MPI_Cart_sub(config.grid_comm, config.coords, &config.col_comm);
	MPI_Comm_rank(config.col_comm, &config.col_rank);
	if(config.world_rank == 63){
		int source;
		int dest;
			printf("%d\n", config.row_rank);
			printf("%d\n", config.col_rank);
			MPI_Cart_shift(config.grid_comm, 0, 1, &source, &dest);
			printf("%d\n", source);
			printf("%d\n", dest);
			MPI_Cart_shift(config.row_comm, 0, 1, &source, &dest);
			printf("%d\n", source);
			printf("%d\n", dest);
			MPI_Cart_shift(config.col_comm, 0, 1, &source, &dest);
			printf("%d\n", source);
			printf("%d\n", dest);
	}
	//MPI_Bcast(1, 1, MPI_INT, )
}

void init_matmul(char *A_file, char *B_file, char *outfile)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config.world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);

	/* Copy output file name to configuration */
	config.outfile = outfile;
	
	/* Get matrix size header */
	if(config.world_rank == 0){	
		//Matrix A
		MPI_File_open(MPI_COMM_SELF, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.A_file);
		MPI_File_read_at(config.A_file, 0, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		//MPI_File_close(&config.A_file);
		
		//Matrix B
		MPI_File_open(MPI_COMM_SELF, B_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.B_file);
        MPI_File_read_at(config.B_file, 0, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
        //MPI_File_close(&config.B_file);

		
		/* Verify dim of A and B matches for matul and both are square*/
		if(!(config.A_dims[0] == config.B_dims[0] && config.A_dims[1] == config.B_dims[1])){
			MPI_Finalize();
			exit(1);
		}
		
	}
	/* Broadcast global matrix sizes */
	MPI_Bcast(config.A_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(config.B_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	config.matrix_size = config.A_dims[0];
	
	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
	config.dim[0] = sqrt(config.world_size);
	config.dim[1] = sqrt(config.world_size);

	/* Create Cart communicator for NxN processes */
	int wrap[2];
	wrap[0] = wrap[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, wrap, 1, &config.grid_comm);

	/* Sub div cart communicator to N row communicator */
	config.coords[0] = 0;
	config.coords[1] = 1;
	MPI_Cart_sub(config.grid_comm, config.coords, &config.row_comm);
	/* Sub div cart communicator to N col communicator */
	config.coords[0] = 1;
	config.coords[1] = 0;
	MPI_Cart_sub(config.grid_comm, config.coords, &config.col_comm);

	MPI_Comm_rank(config.row_comm, &config.row_rank);
	MPI_Comm_size(config.row_comm, &config.row_size);
	
	MPI_Comm_rank(config.col_comm, &config.col_rank);
	MPI_Comm_size(config.col_comm, &config.col_size);
	
	if(config.world_rank == 0 || config.world_rank == 1 || config.world_rank == 2){
		printf("%d,%d,%d\n", config.row_rank, config.col_rank, config.world_rank);
	}

	/* Setup sizes of full matrices */

	/* Setup sizes of local matrix tiles */
	config.local_size = config.matrix_size / sqrt(config.world_size);
	config.local_dims[0] = config.local_size;
	config.local_dims[1] = config.local_size;

	/* Create subarray datatype for local matrix tile */

	int startOffset[2] = {config.local_size * config.row_rank,config.local_size*config.col_rank};

	MPI_Datatype toTile;
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, startOffset, MPI_ORDER_C, MPI_DOUBLE, &toTile);
	MPI_Type_commit(&toTile);

	
	/* Create data array to load actual block matrix data */
	//double matrixData[10*10];
	config.A = malloc(sizeof(double) * (config.local_size * config.local_size));
	config.B = malloc(sizeof(double) * (config.local_size * config.local_size));
	//config.C = malloc(sizeof(double) * (config.A_dims[0] * config.A_dims[1]));


	/* Set fileview of process to respective matrix block */

	MPI_Offset offset = 2 * sizeof(int);
	MPI_File_set_view(config.A_file, offset , toTile, MPI_DOUBLE , "native", MPI_INFO_NULL);
	MPI_File_set_view(config.B_file, offset , toTile, MPI_DOUBLE , "native", MPI_INFO_NULL);
		
	/* Collective read blocks from files */
	
	MPI_File_read_all(config.A_file, config.A, config.local_size*config.local_size,MPI_DOUBLE,  MPI_STATUS_IGNORE);
	MPI_File_read_all(config.B_file, config.B, config.local_size*config.local_size,MPI_DOUBLE,  MPI_STATUS_IGNORE);


	/* Close data source files */
	MPI_File_close(&config.A_file);
	MPI_File_close(&config.B_file);
	MPI_Type_free(&toTile);
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

	/* Compute source and target for vertical shift of B blocks */
	int source, dest;
	MPI_Cart_shift(config.row_comm, 0, 1, &source, &dest);
	int i;
	int root;
	config.A_tmp = malloc(sizeof(config.local_size));
	for (i = 0; i < config.dim[0]; i++) {
		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */
		root = (config.row_rank + i) % config.dim[0];
		if(root == config.col_rank){
			MPI_Bcast(config.A, 1, config.block, root, config.col_comm);

		} else{
			MPI_Bcast(config.A_tmp, 1, config.block, root, config.row_comm);
		}
		MPI_Sendrecv_replace(config.B, 1, config.block, dest, 0, source, 0, config.col_comm, MPI_STATUS_IGNORE);
		/*if(i == config.row_coll && i == config.row_rank){

		}
		MPI_Bcast(config.A, 2, MPI_INT, broad, config.row_comm);*/
		/* dgemm with blocks */
		
		/* Shfting block B upwards and receive from process below */
		//MPI_Cart_shift(config.col_comm, 0, 1, &source, &dest);
	}
}
