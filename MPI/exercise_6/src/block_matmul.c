#include "block_matmul.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

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

	/* Copy output file name to configuration */
	config.outfile = outfile;
	
	/* Get matrix size header */
	MPI_File_open(MPI_COMM_SELF, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.A_file);
	MPI_File_open(MPI_COMM_SELF, B_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.B_file);
	if(config.world_rank == 0){
		
		//Matrix A
		//MPI_File_open(MPI_COMM_SELF, A_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.A_file);
		MPI_File_read_at(config.A_file, 0, config.A_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		//MPI_File_close(&config.A_file);
		
		//Matrix B
        MPI_File_read_at(config.B_file, 0, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
        //MPI_File_close(&config.B_file);

		
		/* Verify dim of A and B matches for matul and both are square*/
		if(!(config.A_dims[0] == config.B_dims[0] && config.A_dims[1] == config.B_dims[1])){
			MPI_Finalize();
			exit(1);
		}

		MPI_Offset tempoffset = 2 * sizeof(int);
		double temp[config.A_dims[0] * config.A_dims[0]];
		MPI_File_read_at(config.A_file, tempoffset, &temp, config.A_dims[0]*config.A_dims[0], MPI_DOUBLE, MPI_STATUS_IGNORE);
		int i,j;
		printf("\n");
		printf("-----FULL MATRIX-----\n");
		for(i = 0; i < config.A_dims[0]; i++){
			for(j = 0; j < config.A_dims[0]; j++){
				printf("%f  ", temp[i*config.A_dims[0]+j]);
			}
			printf("\n");
		}
		printf("-----FULL MATRIX-----\n");
		printf("\n");
		
	}
	/* Broadcast global matrix sizes */
	MPI_Bcast(config.A_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(config.B_dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	config.matrix_size = config.A_dims[0];
	//config.C_dims = config.A_dims;
	config.C_dims[0] = config.A_dims[0];
	config.C_dims[1] = config.A_dims[1];
	
	/* Set dim of tiles relative to the number of processes as NxN where N=sqrt(world_size) */
	config.dim[0] = sqrt(config.world_size);
	config.dim[1] = sqrt(config.world_size);

	/* Create Cart communicator for NxN processes */
	int wrap[2];
	wrap[0] = wrap[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, config.dim, wrap, 1, &config.grid_comm);

	/* Sub div cart communicator to N row communicator */
	//config.coords[0] = 0;//original
	//config.coords[1] = 1;//original
	//MPI_Cart_sub(config.grid_comm, config.coords, &config.row_comm);
	int klax[2] = {0,1};
	MPI_Cart_sub(config.grid_comm, klax, &config.row_comm);
	/* Sub div cart communicator to N col communicator */
	//config.coords[0] = 1;//original
	//config.coords[1] = 0;//original
	//MPI_Cart_sub(config.grid_comm, config.coords, &config.col_comm);
	int krax[2] = {1,0};
	MPI_Cart_sub(config.grid_comm, krax, &config.col_comm);

	MPI_Comm_rank(config.row_comm, &config.row_rank);
	MPI_Comm_size(config.row_comm, &config.row_size);
	
	MPI_Comm_rank(config.col_comm, &config.col_rank);
	MPI_Comm_size(config.col_comm, &config.col_size);

	/* Setup sizes of full matrices */

	/* Setup sizes of local matrix tiles */
	config.local_size = config.matrix_size / sqrt(config.world_size);
	//printf("local size: %d\n", config.local_size);
	config.local_dims[0] = config.local_size;
	config.local_dims[1] = config.local_size;

	/* Create subarray datatype for local matrix tile */

	int startOffset[2] = {config.local_size*config.col_rank, config.local_size * config.row_rank};
	//startOffset[0] = 5;
	//startOffset[1] = 5;
	//printf("kallestropp %d %d %d %d\n", config.A_dims[0], config.A_dims[1], config.local_dims[0], config.local_dims[1]);
	MPI_Type_create_subarray(2, config.A_dims, config.local_dims, startOffset, MPI_ORDER_C, MPI_DOUBLE, &config.block);
	MPI_Type_commit(&config.block);

	
	/* Create data array to load actual block matrix data */
	config.A = (double *)malloc(sizeof(double) * (config.local_size * config.local_size));
	config.A_tmp = (double *)malloc(sizeof(double) * (config.local_size * config.local_size));
	config.B = (double *)malloc(sizeof(double) * (config.local_size * config.local_size));
	config.C = (double *)malloc(sizeof(double) * (config.local_size * config.local_size));

	int i;
	for(i = 0; i < config.local_size * config.local_size; i++){
		config.A[i] = 0.0;
		config.A_tmp[i] = 0.0;
		config.B[i] = 0.0;
		config.C[i] = 0.0;
	}


	/* Set fileview of process to respective matrix block */

	MPI_Offset offset = 2 * sizeof(int);

	//offset += ((config.world_rank % config.col_size) * config.local_size + config.col_size * config.col_rank * config.local_size*config.local_size) * sizeof(double);
	//MPI_File_set_view(config.A_file, offset , config.block, MPI_DOUBLE , "native", MPI_INFO_NULL);
	MPI_File_set_view(config.A_file, offset , MPI_DOUBLE, config.block , "native", MPI_INFO_NULL);
	MPI_File_set_view(config.B_file, offset , config.block, MPI_DOUBLE , "native", MPI_INFO_NULL);
		
	/* Collective read blocks from files */
	int test = config.local_size*config.local_size;
	MPI_File_read_all(config.A_file, config.A, test , MPI_DOUBLE,  MPI_STATUS_IGNORE);
	MPI_File_read_all(config.B_file, config.B, config.local_size*config.local_size,MPI_DOUBLE,  MPI_STATUS_IGNORE);


	/* Close data source files */
	MPI_File_close(&config.A_file);
	MPI_File_close(&config.B_file);
	//MPI_Type_free(&toTile);
	//time_t t;
	//srand((unsigned) time(&t));
	//sleep((int)(rand()*3));
	sleep(config.world_rank);
	int j;
	printf("\n");
	printf("[%d]\n]", config.world_rank);
	printf("-----Tile-----\n");
	for(i = 0; i < config.local_size; i++){
		for(j = 0; j < config.local_size; j++){
			printf("%f  ", config.A[i*config.local_size+j]);
		}
		printf("\n");
	}
	printf("-----Tile-----\n");
	printf("\n");
}

void cleanup_matmul()
{

	printf("struts\n");
	printf("pasta: %d\n", config.C_dims[0]);

	/* Rank zero writes header specifying dim of result matrix C */
	MPI_File_open(MPI_COMM_SELF, config.outfile, MPI_MODE_RDWR, MPI_INFO_NULL, &config.C_file);
	if(config.world_rank == 0){
		MPI_File_write_at(config.C_file, 0, config.C_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_read_at(config.C_file, 0, config.C_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
	}
	printf("pirre: %d\n", config.C_dims[0]);
	/* Set fileview of process to respective matrix block with header offset */
	MPI_Offset offset = 2 * sizeof(int);
	MPI_File_set_view(config.C_file, offset , config.block, MPI_DOUBLE , "native", MPI_INFO_NULL);

	/* Collective write and close file */
	MPI_File_write_all(config.C_file, config.C, config.local_size*config.local_size, MPI_DOUBLE,  MPI_STATUS_IGNORE);
	
	if(config.world_rank == 0){
		int jump = sizeof(double) * 8;
		double temp;
		MPI_File_read_at(config.C_file, jump, &temp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
		MPI_File_read_at(config.C_file, 0, config.C_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		printf("temp %f\n", temp);
		printf("c dim 1 %d\n", config.C_dims[1]);
		MPI_File_close(&config.C_file);

		MPI_File testFile;
		MPI_File_open(MPI_COMM_SELF, config.outfile, MPI_MODE_RDWR, MPI_INFO_NULL, &testFile);
		MPI_File_read_at(testFile, jump, &temp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
		MPI_File_read_at(testFile, 0, config.C_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		printf("andra temp %f\n", temp);
		printf("andra c dim 1 %d\n", config.C_dims[1]);
		MPI_File_close(&testFile);
	}

	MPI_File_close(&config.C_file);
	/* Cleanup */
}


void compute_fox()
{	

	/* Compute source and target for vertical shift of B blocks */
	int source, dest;
	
	int tileSize = config.local_size * config.local_size;
	MPI_Cart_shift(config.col_comm, 0, 1, &source, &dest);
	int rootX = config.col_rank;
	int i;
	for (i = 0; i < config.dim[0]; i++) {

		int rowID;
		int colID;
		int inRow;

		MPI_Comm_rank(config.row_comm, &rowID);
		MPI_Comm_rank(config.col_comm, &colID);
		MPI_Comm_size(config.row_comm, &inRow);
		
		double **AMul;
		if(rootX == config.row_rank){
			AMul = &config.A;
			printf("upper: [%d]   RowID:%d   ColID:%d   RootX:%d\n", config.world_rank, rowID, colID, rootX);
		}else{
			AMul = &config.A_tmp;
			printf("lower: [%d]   RowID:%d   ColID:%d   RootX:%d\n", config.world_rank, rowID, colID, rootX);
		}
		
		MPI_Bcast(*AMul, tileSize, MPI_DOUBLE, rootX, config.row_comm);

		//AMul = &config.A;
		
		int a,b,c;
		for(a=0;a<config.local_size;a++){
			for(b=0;b<config.local_size;b++){
				int indexC = a * config.local_size + b;
				//config.C[indexC] = 0;
				for(c=0;c<config.local_size;c++){
					int indexA = a * config.local_size + c;
					int indexB = c * config.local_size + b;
					//config.C[indexC]+=config.A[indexA]*config.B[indexB];
					config.C[indexC]+=(*AMul)[indexA]*config.B[indexB];

					double val1 = config.A[indexA];
					double val2 = config.B[indexB];
					//if(i > 0){
					//if(rootX == config.row_rank){
						//printf("[%d] val1: %f   val2:%f\n", config.world_rank, val1, val2);
					//}
					
				}
				//printf("snopp %f\n", config.C[indexC]);
			}
		}
		

		//printf("kraxflax\n");
		//sleep(5);

		//MPI_Cart_shift(config.col_comm, 0, 1, &source, &dest);
		//MPI_Sendrecv_replace(config.B, tileSize, MPI_DOUBLE, dest, config.col_rank, source, source, config.col_comm, MPI_STATUS_IGNORE);
		MPI_Cart_shift(config.col_comm, 0, 1, &source, &dest);
		if(config.world_rank == 0){
			printf("[%d]   source:%d    dest:%d \n", config.world_rank, source, dest);
		}
		MPI_Sendrecv_replace(config.B, tileSize, MPI_DOUBLE, dest, 0, source, 0, config.col_comm, MPI_STATUS_IGNORE);

		rootX = (rootX+1)%config.row_size;

	}
	
	for(i = 0; i < tileSize; i++){
		if(config.C[i] != 0){
			printf("[%d]kalle %f\n",config.world_rank, config.C[i]);
		}
		
	}
}
