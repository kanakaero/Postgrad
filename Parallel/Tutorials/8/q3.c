// Code to demonstrate the use of #pragma acc serial, num_gangs, num_workers, vector_length, private, firstprivate and if clauses in OpenACC.

#include <stdio.h>
#include <stdlib.h>

int main()
{
    int N, i;

    N = 1;

#pragma acc parallel num_gangs(4) num_workers(2) vector_length(8) private(i) firstprivate(N)
    {
        printf("Value of N is %d \n", N);
#pragma acc loop gang
        for (i = 0; i < N; i++)
            printf("Hello from gang! \n");
#pragma acc loop worker
        for (i = 0; i < N; i++)
            printf("Hello from worker! \n");
#pragma acc loop vector
        for (i = 0; i < N; i++)
            printf("Hello from vector! \n");
    }

#pragma acc parallel if (N > 0)
    {
        printf("This parallel region is valid! \n");
    }

    return 0;
}