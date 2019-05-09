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

	//MPI_Request request[num_ranks-1];
	//MPI_Status status[num_ranks-1];

	double *temps = NULL;
	if(world_rank == 0){
		temps = malloc(sizeof(double) * num_ranks);
	}

	MPI_Gather(local_count, 1, MPI_INT, temps, num_ranks, MPI_INT, 0, MPI_COMM_WORLD);

	if (world_rank == 0) {
		double count = 0;
		for(i = 0; i < num_ranks; i++){
			count += temps[i];
		}
		
		double P = (double)count / (double)flip;
		double pi = 4 * P;

		*answer = pi;
	}
	
}
