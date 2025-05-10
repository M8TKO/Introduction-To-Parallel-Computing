#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main(int argc, char **argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_procs = 16;
    integer n = 1000;

    int small = n / num_procs;
    int large = small + 1;
    int num_large = n % num_procs;
    int first, last;
    integer block = (rank < num_large ? large : small);

    double *A = calloc(n * n, sizeof(double));
    double *B = calloc(n * n, sizeof(double));
    double *C = calloc(n * n, sizeof(double));
    double *A_scatter = malloc(n * block * sizeof(double));
    double *A_temp = malloc(n * large * sizeof(double));
    double *B_scatter = malloc(n * block * sizeof(double));
    double *C_local = malloc(n * block * sizeof(double));
    int *sendcounts = malloc(num_procs * sizeof(int));
    int *displs = malloc(num_procs * sizeof(int));

    // Set up sendcounts and displacements for scatter/gather
    for (int i = 0; i < num_procs; i++) {
        sendcounts[i] = n * (i < num_large ? large : small);
        displs[i] = n * (i < num_large ? i * large : num_large * large + (i - num_large) * small);
    }

    // Initialize matrices A and B with random values on rank 0
    if (rank == 0) {
        integer idist = 3, iseed[4] = {3, 54, 244, 7};
        integer N = n * n;
        dlarnv_(&idist, iseed, &N, A);
        dlarnv_(&idist, iseed, &N, B);
    }

    // Scatter matrix A and B
    MPI_Scatterv(A, sendcounts, displs, MPI_DOUBLE, A_scatter, n * block, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, sendcounts, displs, MPI_DOUBLE, B_scatter, n * block, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Initialize C_local to zero
    doublereal alpha = 0;
    dlaset_("A", &n, &block, &alpha, &alpha, C_local, &n);

    // Blocked matrix multiplication
    alpha = 1;
    for (int i = 0; i < num_procs; i++) {
        if (rank == i) {
            dlacpy_("A", &n, &block, A_scatter, &n, A_temp, &n);
        }

        MPI_Bcast(A_temp, n * large, MPI_DOUBLE, i, MPI_COMM_WORLD);
        integer k = sendcounts[i] / n;
        dgemm_("N", "N", &n, &block, &k, &alpha, A_temp, &n, &B_scatter[displs[i] / n], &n, &alpha, C_local, &n);
    }

    // Gather the final result matrix C
    MPI_Gatherv(C_local, n * block, MPI_DOUBLE, C, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Verify result using sequential multiplication (rank 0 only)
    if (rank == 0) {
        doublereal *Reference = malloc(n * n * sizeof(doublereal));
        doublereal beta = 0;
        alpha = 1;
        dgemm_("N", "N", &n, &n, &n, &alpha, A, &n, B, &n, &beta, Reference, &n);

        int i;
        for (i = 0; i < n * n; i++) {
            if (abs(C[i] - Reference[i]) > 1e-12)
                break;
        }

        printf("Result is %svalid!\n", (i == n * n) ? "" : "not ");
        free(Reference);
    }

    free(A);
    free(B);
    free(C);
    free(A_scatter);
    free(A_temp);
    free(B_scatter);
    free(C_local);
    free(sendcounts);
    free(displs);

    MPI_Finalize();
    return 0;
}
