#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TYPE float
#define SMALLVALUE 0.001

double get_time()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// Allocate 2D array
TYPE **alloc_matrix(int N)
{
    TYPE **mat = (TYPE **)malloc(N * sizeof(TYPE *));
    for (int i = 0; i < N; i++)
        mat[i] = (TYPE *)malloc(N * sizeof(TYPE));
    return mat;
}

void free_matrix(TYPE **mat, int N)
{
    for (int i = 0; i < N; i++)
        free(mat[i]);
    free(mat);
}

// Initialize symmetric positive definite matrix
void init(TYPE **mat, int N)
{
    for (int ii = 0; ii < N; ++ii)
        for (int jj = 0; jj < ii; ++jj)
        {
            mat[ii][jj] = (ii + jj) / (float)(N * N);
            mat[jj][ii] = mat[ii][jj];
        }

    for (int ii = 0; ii < N; ++ii)
        mat[ii][ii] = 1.0;
}

// Serial Cholesky
void cholesky_serial(TYPE **a, int N)
{
    for (int ii = 0; ii < N; ++ii)
    {
        for (int jj = 0; jj < ii; ++jj)
        {
            for (int kk = 0; kk < jj; ++kk)
                a[ii][jj] -= a[ii][kk] * a[jj][kk];

            a[ii][jj] /= (a[jj][jj] > SMALLVALUE ? a[jj][jj] : 1);
        }

        for (int kk = 0; kk < ii; ++kk)
            a[ii][ii] -= a[ii][kk] * a[ii][kk];

        a[ii][ii] = sqrt(a[ii][ii]);
    }
}

// OpenACC Cholesky
void cholesky_acc(TYPE **a, int N)
{
#pragma acc data copy(a[0:N][0:N])
    {
        for (int ii = 0; ii < N; ++ii)
        {
            for (int jj = 0; jj < ii; ++jj)
            {
                TYPE sum = 0.0;

#pragma acc parallel loop reduction(+ : sum)
                for (int kk = 0; kk < jj; ++kk)
                    sum += a[ii][kk] * a[jj][kk];

                a[ii][jj] = (a[ii][jj] - sum) / a[jj][jj];

#pragma acc update device(a[ii][jj])
            }

            TYPE sum = 0.0;

#pragma acc parallel loop reduction(+ : sum)
            for (int kk = 0; kk < ii; ++kk)
                sum += a[ii][kk] * a[ii][kk];

            a[ii][ii] = sqrt(a[ii][ii] - sum);

#pragma acc update device(a[ii][ii])
        }
    }
}

// Compare two lower triangular matrices
double compute_error(TYPE **a, TYPE **b, int N)
{
    double err = 0.0;
    double norm = 0.0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j <= i; j++)
        {
            double diff = a[i][j] - b[i][j];
            err += diff * diff;
            norm += a[i][j] * a[i][j];
        }

    return sqrt(err) / sqrt(norm);
}

// Validate A ≈ L L^T
double validate_cholesky(TYPE **L, TYPE **A_orig, int N)
{
    TYPE **A_rec = alloc_matrix(N);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A_rec[i][j] = 0.0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j <= i; j++)
        {
            for (int k = 0; k <= j; k++)
                A_rec[i][j] += L[i][k] * L[j][k];
            A_rec[j][i] = A_rec[i][j];
        }

    double err = 0.0, norm = 0.0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            double diff = A_rec[i][j] - A_orig[i][j];
            err += diff * diff;
            norm += A_orig[i][j] * A_orig[i][j];
        }

    free_matrix(A_rec, N);

    double rel_error = sqrt(err) / sqrt(norm);
    return rel_error;
}

int main()
{
    FILE *fp = fopen("cholesky_results.dat", "w");
    if (!fp)
    {
        perror("File open failed");
        return 1;
    }

    fprintf(fp, "# N SerialTime OpenACCTime RelError ReconError\n");

    int sizes[] = {10, 100, 1000};
    int num_sizes = 3;

    for (int s = 0; s < num_sizes; s++)
    {
        int N = sizes[s];

        TYPE **a = alloc_matrix(N);
        TYPE **b = alloc_matrix(N);
        TYPE **A_orig = alloc_matrix(N);

        init(a, N);

        // copy
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                b[i][j] = a[i][j];
                A_orig[i][j] = a[i][j];
            }

        double t1 = get_time();
        cholesky_serial(a, N);
        double t2 = get_time();

        double t3 = get_time();
        cholesky_acc(b, N);
        double t4 = get_time();

        double error = compute_error(a, b, N);
        double recon_error = validate_cholesky(b, A_orig, N);

        printf("N = %d\n", N);
        printf("Serial time   : %.6f s\n", t2 - t1);
        printf("OpenACC time  : %.6f s\n", t4 - t3);
        printf("Relative error: %.6e\n", error);
        printf("Recon error   : %.6e\n\n", recon_error);

        fprintf(fp, "%d %.6f %.6f %.6e %.6e\n",
                N, t2 - t1, t4 - t3, error, recon_error);

        free_matrix(a, N);
        free_matrix(b, N);
        free_matrix(A_orig, N);
    }

    fclose(fp);
    return 0;
}