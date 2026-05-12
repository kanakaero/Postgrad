#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <curand.h>

#define N 100
#define PROBZERO 33

int main()
{
    int a[N][N];
    int x[N] = {0}, y[N] = {0};

    curandGenerator_t gen;
    curandCreateGeneratorHost(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, time(NULL));

    // Generate random numbers in bulk
    int *rand_nums = (int *)malloc(N * N * sizeof(int));
    curandGenerate(gen, (unsigned int *)rand_nums, N * N);


# pragma acc parallel loop collapse(2) copyin(rand_nums[ : N * N]) copyout(a[ : N][ : N])
    for (int ii = 0; ii < N; ++ii)
    {
        for (int jj = 0; jj < N; ++jj)
        {
            int idx = ii * N + jj;
            int r = rand_nums[idx] % N;
            if ((rand_nums[idx] % 100) < PROBZERO)
                r = 0;

            a[ii][jj] = r;
        }
    }

#pragma acc parallel loop collapse(2) reduction(+ : x[ : N]) reduction(+ : y[ : N])
    for (int ii = 0; ii < N; ++ii)
    {
        for (int jj = 0; jj < N; ++jj)
        {
            x[ii] += a[ii][jj];
            if (a[ii][jj] > 0)
                y[ii]++;
        }
    }

    for (int ii = 0; ii < N; ++ii)
        printf("%.0f ", x[ii] * 100.0 / (y[ii] * N));

    printf("\n");

    free(rand_nums);
    curandDestroyGenerator(gen);
    return 0;
}