#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int MYPROC = rank + 1;

    int n = 100;
    int num_procs = 15;
    int proc_n = num_procs + 1;
    int num_leaves = proc_n / 2;
    int m = num_leaves;
    int k = log2(num_leaves * 2);

    int *proc_ranks;
    MPI_Group world_group, leaves_group;
    MPI_Comm leaves_comm;

    proc_ranks = malloc(num_leaves * sizeof(int));
    for (int i = 0; i < num_leaves; i++)
        proc_ranks[i] = num_leaves - 1 + i;

    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, num_leaves, proc_ranks, &leaves_group);
    MPI_Comm_create(MPI_COMM_WORLD, leaves_group, &leaves_comm);

    double *x, *y, *local_x, *local_y, S, L, R;

    // Leaf nodes perform local prefix sum
    if (rank >= num_leaves - 1) {
        int block = (n % num_leaves ? n / num_leaves + 1 : n / num_leaves);
        int block_dim = num_leaves * block;
        local_x = malloc(block * sizeof(double));
        local_y = malloc(block * sizeof(double));

        if (rank == num_leaves - 1) {
            x = malloc(block_dim * sizeof(double));
            y = malloc(block_dim * sizeof(double));
            for (int i = 0; i < n; i++)
                x[i] = 1;
            for (int i = n; i < block_dim; i++)
                x[i] = 0;
        }

        MPI_Scatter(x, block, MPI_DOUBLE, local_x, block, MPI_DOUBLE, 0, leaves_comm);

        local_y[0] = local_x[0];
        for (int i = 1; i < block; i++)
            local_y[i] = local_y[i - 1] + local_x[i];

        S = local_y[block - 1];
        MPI_Send(&S, 1, MPI_DOUBLE, MYPROC / 2 - 1, 0, MPI_COMM_WORLD);

        // Wait for downward step
        if (rank > num_leaves - 1) {
            MPI_Recv(&R, 1, MPI_DOUBLE, MYPROC / 2 - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < block; i++)
                local_y[i] += R;
        }

        // Leftmost leaf gathers and prints result
        MPI_Gather(local_y, block, MPI_DOUBLE, y, block, MPI_DOUBLE, 0, leaves_comm);

        if (rank == num_leaves - 1) {
            for (int i = 0; i < n; i++)
                printf("%g\n", y[i]);
            free(x);
            free(y);
        }

        free(local_x);
        free(local_y);
    }

    // Internal nodes: Upward phase
    for (int j = 2; j <= k - 1; j++) {
        m /= 2;
        if (MYPROC >= m && MYPROC <= 2 * m - 1) {
            MPI_Recv(&L, 1, MPI_DOUBLE, 2 * MYPROC - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&R, 1, MPI_DOUBLE, 2 * MYPROC, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            S = L + R;
            MPI_Send(&S, 1, MPI_DOUBLE, MYPROC / 2 - 1, 0, MPI_COMM_WORLD);
        }
    }

    // Root: transition from up to down phase
    m = 1;
    if (MYPROC == 1) {
        MPI_Recv(&L, 1, MPI_DOUBLE, 2 * MYPROC - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&R, 1, MPI_DOUBLE, 2 * MYPROC, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        S = L + R;
        MPI_Send(&L, 1, MPI_DOUBLE, 2 * MYPROC, 0, MPI_COMM_WORLD);
    }

    // Downward phase
    for (int j = 1; j <= k - 2; j++) {
        m *= 2;
        if (MYPROC == m) {
            MPI_Send(&L, 1, MPI_DOUBLE, 2 * MYPROC, 0, MPI_COMM_WORLD);
        } else if (MYPROC >= m + 1 && MYPROC <= 2 * m - 1) {
            MPI_Recv(&R, 1, MPI_DOUBLE, MYPROC / 2 - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&R, 1, MPI_DOUBLE, 2 * MYPROC - 1, 0, MPI_COMM_WORLD);
            S = L + R;
            MPI_Send(&S, 1, MPI_DOUBLE, 2 * MYPROC, 0, MPI_COMM_WORLD);
        }
    }

    free(proc_ranks);
    MPI_Finalize();
    return 0;
}
