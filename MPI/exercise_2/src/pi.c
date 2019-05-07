#include "pi.h"
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

	double radius = 1.0;
	int i;
	for(i = 0; i < flip / (double)num_ranks; i++){
		double newCoordX = (double)rand()/RAND_MAX;
		double newCoordY = (double)rand()/RAND_MAX;

		double distance = sqrt(newCoordX*newCoordX+newCoordY*newCoordY);

		if(distance <= radius){
			*local_count++;
		}
	}
	
	if (world_rank == 0) {
		int count = *local_count;
		int temp;
		for(i = 1; i < num_ranks; i++){
			
    		MPI_Recv(&temp, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			count += temp;
		}
		
		double P = (double)count / (double)flip;
		double pi = 4 * P;

		*answer = pi;

		printf("pi: %f\n", pi);
	}else{
		MPI_Send(local_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	
}
