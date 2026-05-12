#include <stdio.h>
#include <stdlib.h>

int main()
{
    int mynumber, i, N, *numbers, flag = 0, iloc=0;

    numbers = (int *)malloc(10 * sizeof(int));
    numbers[0] = 1;
    numbers[1] = 2;
    numbers[2] = 3;
    numbers[3] = 4;
    numbers[4] = 5;
    numbers[5] = 6;
    numbers[6] = 7;
    numbers[7] = 8;
    numbers[8] = 9;
    numbers[9] = 10;

    N = 10;
    printf("Enter the number to be searched \n");
    scanf("%d", &mynumber);

#pragma acc parallel num_gangs(5) firstprivate(N) copyin(numbers[0:N]) copy(flag) copy(iloc)
    {
        #pragma acc loop gang reduction(max:flag)
        for (i = 0; i < N; i++)
        {
            if (numbers[i] == mynumber)
            {
                flag = 1;
                #pragma acc atomic write
                iloc = i;
            }
        }
    }

    if (flag == 1)
        printf("Number found at location %d \n", iloc);
    else
        printf("Number not found \n");
}