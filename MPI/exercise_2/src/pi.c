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
	srand(world_rank);
	double radius = 1.0;
	int i;
	for(i = 0; i < (int)((double)flip / (double)num_ranks); i++){
		double newCoordX = (double)rand()/RAND_MAX;
		double newCoordY = (double)rand()/RAND_MAX;
		//printf("%f\n", newCoordX);
		double distance = sqrt(newCoordX*newCoordX+newCoordY*newCoordY);
		
		if(distance <= radius){
			//printf("%d   ", *local_count);
			(*local_count)++;
			//printf("%d\n", *local_count);
		}
	}

	MPI_Request request[num_ranks-1];
	MPI_Status status[num_ranks-1];
	
	if (world_rank == 0) {
		int count = 0;
		int temp[num_ranks];
		temp[0] = *local_count;



		for(i = 1; i < num_ranks; i++){
			
    		//MPI_Irecv(&(temp[i]), 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Irecv(&(temp[i]), 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		MPI_Waitall(num_ranks-1, request, status);

		for(i = 0; i < num_ranks; i++){
			count += temp[i];
		}
		
		double P = (double)count / (double)flip;
		double pi = 4 * P;

		*answer = pi;
	}else{
		//MPI_Send(local_count, 1, MPI_INT, 0, world_rank, MPI_COMM_WORLD);
		//&data, count, MPI_INT, source, tag, MPI_COMM_WORLD, &request[1]
		MPI_Send(local_count, 1, MPI_INT, 0, world_rank, MPI_COMM_WORLD, &(request[world_rank]));
	}
	
}
