#include "block_matmul.h"
#include <unistd.h>


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
        MPI_Comm_rank(MPI_COMM_WORLD, &config.world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	if(config.world_rank == 0){
		printf("kalle");
		MPI_File_open(MPI_COMM_SELF, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.A_file);
                MPI_File_read_at(config.A_file, 0, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		printf("%d  %d\n", config.A_dims[0], config.A_dims[1]);
		
		double krabba[25];
		//int offset = sizeof(int)*2;
		MPI_Offset offset = 2 * sizeof(int);
		MPI_File_read_at(config.A_file, offset, krabba, 25, MPI_DOUBLE, MPI_STATUS_IGNORE);
		//printf("%f  %f\n", krabba[0], krabba[1]);
		printf("BigArr:\n");
		int i,j;
		/*
		for(i = 0; i < 5; i++){
                	for(j = 0; j < 5; j++){
				printf("%f ", krabba[i*5+j]);
                	}
			printf("\n");
		}*/
		for(i = 0; i < 25; i++){
			printf("%f ", krabba[i]);
		}
		printf("\n");
		//int starts[2] = {3,1};
		int starts[1] = {7};
		int littleSize[1] = {10};
		int bigSize[1] = {25};
		
		MPI_Datatype mysubarray;
		MPI_Type_create_subarray(1, bigSize, littleSize, starts, MPI_ORDER_C, MPI_DOUBLE, &mysubarray);
		MPI_Type_commit(&mysubarray);
		
		MPI_Request pelle;
		MPI_Isend(krabba, 1, mysubarray, config.world_rank, config.world_rank, MPI_COMM_WORLD, &pelle);
		double littleKrabba[10];
		//MPI_Wait(&pelle, MPI_STATUS_IGNORE);
		MPI_Recv(littleKrabba, 10, mysubarray, config.world_rank, config.world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
		printf("LittleArr:\n");
		for(i = 0; i < 10; i++){
                        printf("%f ", littleKrabba[i]);
                }
                printf("\n");

		MPI_Type_free(&mysubarray);
                MPI_File_close(&config.A_file);
		
	}
}

void init_matmul_snurresprett(char *A_file, char *B_file, char *outfile)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config.world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &config.world_size);
	MPI_Comm lastBarrier;
	MPI_Comm_split(MPI_COMM_WORLD, 0, config.world_rank, &lastBarrier);

	//MPI_Bcast(&lastBarrier, 1, MPI_Comm, 0, MPI_COMM_WORLD);
	//MPI_Bcast(&lastBarrier, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//MPI_Request request;
	if(config.world_rank == 0){
		printf("World size: %d\n", config.world_size);
		printf("thread: %d\n", config.world_rank);
	}else{
		//MPI_Wait(&request, MPI_STATUS_IGNORE);
		MPI_Recv(NULL, 0, MPI_INT, config.world_rank-1, config.world_rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("thread: %d\n", config.world_rank);
	}
	if(config.world_rank != config.world_size-1){
        	MPI_Send(NULL, 0, MPI_INT, config.world_rank+1, config.world_rank, MPI_COMM_WORLD);
        }
	if(config.world_rank == 10){
		usleep(1000*1000);
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(lastBarrier);
	printf("wtf");
}

void init_matmul_copy(char *A_file, char *B_file, char *outfile)
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
		MPI_File_close(&config.A_file);
		
		//Matrix B
		MPI_File_open(MPI_COMM_SELF, B_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.B_file);
        MPI_File_read_at(config.B_file, 0, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&config.B_file);

		
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

	/* Setup sizes of full matrices */
	config.A = malloc(sizeof(double) * (config.A_dims[0] * config.A_dims[1]));
	config.B = malloc(sizeof(double) * (config.B_dims[0] * config.B_dims[1]));
	config.C = malloc(sizeof(double) * (config.A_dims[0] * config.A_dims[1]));

	/* Setup sizes of local matrix tiles */
	config.local_size = config.matrix_size / sqrt(config.world_size);
	config.local_dims[0] = config.local_size;
	config.local_dims[1] = config.local_size;

	/* Create subarray datatype for local matrix tile */
	config.A_tmp = malloc(sizeof(double) * (config.local_dims[0] * config.local_dims[1]));
	/* Create data array to load actual block matrix data */
	double dataTmp[config.local_dims[0] * config.local_dims[1]];
	/* Set fileview of process to respective matrix block */
	MPI_Offset offset = 2 * sizeof(int);
	int count = config.A_dims[0] * config.A_dims[1]-1;
	MPI_File_open(MPI_COMM_SELF, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.A_file);
	MPI_File_read_at(config.A_file, offset, &dataTmp, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&config.A_file);
	/*MPI_Datatype arraytype;
	MPI_Type_contiguous(3, MPI_DOUBLE, &arraytpe);
	MPI_File_set_view(config.A_file, offset, MPI_DOUBLE, arraytype, "native", MPI_INFO_NULL);
	*/
	printf("%f\n", dataTmp[1]); 
	printf("%f\n", dataTmp[2]); 
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
