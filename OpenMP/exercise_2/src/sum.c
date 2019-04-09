#include "sum.h"

void omp_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum = 0;
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < length(sum_ret); i++){
            sum += sum_ret[i];
        }
    }

}

void omp_critical_sum(double *sum_ret)
{

}

void omp_atomic_sum(double *sum_ret)
{

}

void omp_local_sum(double *sum_ret)
{

}

void omp_padded_sum(double *sum_ret)
{

}

void omp_private_sum(double *sum_ret)
{

}

void omp_reduction_sum(double *sum_ret)
{

}
