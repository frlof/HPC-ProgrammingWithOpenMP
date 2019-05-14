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
		printf("%d\n", config.A_dims[0]);
		printf("%d\n", config.A_dims[1]);
		//Matrix B
		MPI_File_open(MPI_COMM_SELF, B_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &config.B_file);
        MPI_File_read_at(config.B_file, 0, config.B_dims, 2, MPI_INT, MPI_STATUS_IGNORE);
		printf("%d\n", config.A_dims[0]);
		printf("%d\n", config.A_dims[1]);
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
	//printf("krabba %d\n", config.local_size);
	config.A = malloc(sizeof(double) * (config.local_size * config.local_size));
	config.A_tmp = malloc(sizeof(double) * (config.local_size * config.local_size));
	config.B = malloc(sizeof(double) * (config.local_size * config.local_size));
	config.C = malloc(sizeof(double) * (config.local_size * config.local_size));


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
	//int rowID;
	//int inRow;

	//MPI_Comm_rank(config.row_comm, &rowID);
	//MPI_Comm_size(config.row_comm, &inRow);
	
	//print("[%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);



	/* Compute source and target for vertical shift of B blocks */
	int source, dest;
	//MPI_Cart_shift(config.col_comm, 0, 1, &source, &dest);
	
	int tileSize = config.local_size * config.local_size;
	double temp;
	if(config.world_rank == 0){
		temp = 1;
	}
	temp = 0;
	MPI_Bcast(&temp, tileSize, MPI_DOUBLE, 0, config.row_comm);
	int rootX; 
	//int rootY = config.row_rank;
	int i;
	for (i = 0; i < config.dim[0]; i++) {
		rootX = config.row_rank;
		int rowID;
		int inRow;

		MPI_Comm_rank(config.row_comm, &rowID);
		MPI_Comm_rank(config.col_comm, &inRow);
		//printf("localSize: %d\n", config.local_size);
		/*if(config.world_rank == 0){
			printf("[%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
		}*/
		

		double **AMul;
		if(rootX == config.col_rank){
			AMul = &config.A;
		}else{
			AMul = &config.A_tmp;
		}
		printf("[%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
		double temp = 10;
		MPI_Bcast(AMul, tileSize, MPI_DOUBLE, rootX, config.row_comm);
		/*if(rootX == config.col_rank){
			printf("upper: [%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
			MPI_Bcast(config.A, tileSize, MPI_DOUBLE, rootX, config.row_comm);
		}else{
			printf("lower: [%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
			MPI_Bcast(config.A_tmp, tileSize, MPI_DOUBLE, rootX, config.row_comm);
		}*/

		if(config.world_rank == 0){
			printf("[%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
		} 

		if(rootX == config.col_rank){
			/*
			int a,b,c;
			for(a=0;a<config.local_size;a++){
				for(b=0;b<config.local_size;b++){
					for(c=0;c<config.local_size;c++){
						config.C[a][b]+=config.A[a][c]*config.B[c][b];
					}
				}
			}*/
			int a,b,c;
			for(a=0;a<config.local_size;a++){
				for(b=0;b<config.local_size;b++){
					for(c=0;c<config.local_size;c++){

						int indexC = a * config.local_size + b;
						int indexA = a * config.local_size + c;
						int indexB = c * config.local_size + b;
						//config.C[a][b]+=config.A_tmp[a][c]*config.B[c][b];
						config.C[indexC]+=config.A[indexA]*config.B[indexB];
					}
				}
			}

		}else{
			int a,b,c;
			for(a=0;a<config.local_size;a++){
				for(b=0;b<config.local_size;b++){
					int indexC = a * config.local_size + b;
					config.C[indexC] = 0;
					for(c=0;c<config.local_size;c++){

						
						int indexA = a * config.local_size + c;
						int indexB = c * config.local_size + b;
						//config.C[a][b]+=config.A_tmp[a][c]*config.B[c][b];
						config.C[indexC]+=config.A_tmp[indexA]*config.B[indexB];
					}
				}
			}
		}
		//MPI_Request pelle;
		//MPI_Isend(config.B, tileSize, MPI_DOUBLE, 0, config.world_rank, MPI_COMM_WORLD, &pelle);

		//MPI_recv(config.B, tileSize, MPI_INT, i, i, MPI_COMM_WORLD, );
	/*	if(config.world_rank == 0){
			printf("[%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
		}*/
		MPI_Cart_shift(config.row_comm, 0, 1, &source, &dest);

		/*if(config.world_rank == 0){
			printf("[%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
		}*/

		//MPI_Sendrecv(config.B, tileSize, MPI_DOUBLE, dest, sendtag, source, recvtag, comm, status);


		MPI_Sendrecv_replace(config.B, tileSize, MPI_DOUBLE, dest, config.col_rank, source, source, config.col_comm, MPI_STATUS_IGNORE);

	/*	if(config.world_rank == 0){
			printf("[%d]   ID:%d   N:%d\n", config.world_rank, rowID, inRow);
		}*/

		rootX = (rootX+1)%config.row_size;

		



		/* Diag + i broadcast block A horizontally and use A_tmp to preserve own local A */
	//	int broad = 0;

		/*if(i == config.row_coll && i == config.row_rank){

		}
		MPI_Bcast(config.A, 2, MPI_INT, broad, config.row_comm);*/
		/* dgemm with blocks */
		
		/* Shfting block B upwards and receive from process below */
		//MPI_Cart_shift(config.col_comm, 0, 1, &source, &dest);
	}
}
