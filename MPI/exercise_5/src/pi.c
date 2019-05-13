#include "pi.h"

void init_pi(int set_seed, char *outfile)
{
	if (filename != NULL) {
		free(filename);
		filename = NULL;
	}

	if (outfile != NULL) {
		filename = (char*)calloc(sizeof(char), strlen(outfile)+1);
		memcpy(filename, outfile, strlen(outfile));
		filename[strlen(outfile)] = 0;
	}
	seed = set_seed;
}

void cleanup_pi()
{
	if (filename != NULL)
		free(filename);
}

void compute_pi(int flip, int *local_count, double *answer)
{
	int world_rank;
	int num_ranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	srand(world_rank * seed);
	int i;
	for(i = 0; i < (int)((double)flip / (double)num_ranks); i++){
		double x = (double)rand() / (double)RAND_MAX;
		double y = (double)rand() / (double)RAND_MAX;
		double distance = sqrt((x*x) + (y*y));
		if(distance <= 1.0){
			(*local_count)++;
		}
	}
	int count;
	MPI_Reduce(local_count, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (world_rank == 0) {
		
		double P = (double)count / (double)flip;
		double pi = 4 * P;

		*answer = pi;
	}
	MPI_File fh;
	MPI_File_open(MPI_COMM_SELF, "results.out", MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_File_write_at(fh, &local_count, 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
}
