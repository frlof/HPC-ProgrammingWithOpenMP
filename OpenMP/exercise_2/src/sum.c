#include "sum.h"

void omp_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum = 0;
    int iterations = sizeof(sum_ret) / sizeof(sum_ret[0]);
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < iterations; i++){
            sum += sum_ret[i];
        }
    }
}

void omp_critical_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum = 0;
    int iterations = sizeof(sum_ret) / sizeof(sum_ret[0]);
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < iterations; i++){
            #pragma omp critical
            sum += sum_ret[i];
        }
    }
}

void omp_atomic_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum = 0;
    int iterations = sizeof(sum_ret) / sizeof(sum_ret[0]);
    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < iterations; i++){
            #pragma omp atomic
            sum += sum_ret[i];
        }
    }
}

void omp_local_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum[32];
    int iterations = sizeof(sum_ret) / sizeof(sum_ret[0]);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        #pragma omp for
        for(int i = 0; i < iterations; i++){
            sum[id] += sum_ret[i];
        }
    }
}

void omp_padded_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum[32][8];
    int iterations = sizeof(sum_ret) / sizeof(sum_ret[0]);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        #pragma omp for
        for(int i = 0; i < iterations; i++){
            sum[id][0] += sum_ret[i];
        }
    }
}

void omp_private_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum = 0;
    int iterations = sizeof(sum_ret) / sizeof(sum_ret[0]);
    #pragma omp parallel
    {
        int sumPriv = 0;
        #pragma omp for
        for(int i = 0; i < iterations; i++){
            sumPriv += sum_ret[i];
        }
        #pragma omp critical
        sum += sumPriv;
    }
}

void omp_reduction_sum(double *sum_ret)
{
    omp_set_num_threads(32);
    int sum = 0;
    int iterations = sizeof(sum_ret) / sizeof(sum_ret[0]);
    #pragma omp parallel
    {
        #pragma omp for reduction(+:sum)
        for(int i = 0; i < iterations; i++){
            sum += sum_ret[i];
        }
    }
}
