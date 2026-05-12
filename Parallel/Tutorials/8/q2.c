#include <stdio.h>
#include <stdlib.h>

int main()
{
    int mynumber, i, N, *numbers;

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

#pragma acc parallel
    {
        for (i = 0; i < N; i++)
            if (numbers[i] == mynumber)
                break;

        if (i == N)
            printf("Number not found \n");
        else
            printf("Number found at location %d \n", i);
    }
}