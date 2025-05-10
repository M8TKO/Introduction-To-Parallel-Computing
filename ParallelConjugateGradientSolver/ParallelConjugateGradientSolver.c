#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

#define TOLERANCE 1e-8
#define NUM_THREADS 16

// LAPACK-related constants
const integer n = 25;
const doublereal one = 1.0, negOne = -1.0, zero = 0.0;
const integer intOne = 1;
integer INFO;

// Synchronization variables
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
int counter = 0;
int limit = NUM_THREADS;
int stopped = 0;

typedef void (*working_function_t)(int, void*);
working_function_t function;
void* argument;

typedef struct {
    int rank;
} thread_data;

thread_data **wd;

void synchronize() {
    pthread_mutex_lock(&mutex);
    counter++;
    if (counter == limit) {
        counter = 0;
        pthread_cond_broadcast(&cond);
    } else {
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
}

void* thread_loop(void* arg) {
    thread_data* td = (thread_data*) arg;
    while (1) {
        synchronize();
        if (stopped) break;
        function(td->rank, argument);
        synchronize();
    }
    return NULL;
}

// --- Parallel function definitions and their structs ---
typedef struct {
    doublereal *a, *b;
    doublereal alpha, beta;
    doublereal *result;
} vec_op_args;

vec_op_args *vec_op_arg;

void vec_op(int rank, void* arg) {
    vec_op_args *data = (vec_op_args*) arg;
    int small = n / NUM_THREADS;
    int large = small + 1;
    int num_large = n % NUM_THREADS;
    int block = (rank < num_large ? large : small);
    int first = (rank < num_large) ? rank * large : num_large * large + (rank - num_large) * small;
    int last = first + block;

    for (int i = first; i < last; i++) {
        data->result[i] = data->alpha * data->a[i] + data->beta * data->b[i];
    }
}

typedef struct {
    doublereal *src, *dest;
} vec_copy_args;

vec_copy_args *vec_copy_arg;

void vec_copy(int rank, void* arg) {
    vec_copy_args *data = (vec_copy_args*) arg;
    int small = n / NUM_THREADS;
    int large = small + 1;
    int num_large = n % NUM_THREADS;
    int block = (rank < num_large ? large : small);
    int first = (rank < num_large) ? rank * large : num_large * large + (rank - num_large) * small;
    int last = first + block;

    if (last > first) {
        dlacpy_("A", &(integer){last - first}, &intOne, data->src + first, &(integer){last - first}, data->dest + first, &(integer){last - first});
    }
}

typedef struct {
    doublereal *a, *b;
    doublereal *partial_sums;
} dot_args;

dot_args *dot_arg;

void dot_product(int rank, void* arg) {
    dot_args *data = (dot_args*) arg;
    int small = n / NUM_THREADS;
    int large = small + 1;
    int num_large = n % NUM_THREADS;
    int block = (rank < num_large ? large : small);
    int first = (rank < num_large) ? rank * large : num_large * large + (rank - num_large) * small;
    int last = first + block;

    doublereal sum = 0.0;
    for (int i = first; i < last; i++) {
        sum += data->a[i] * data->b[i];
    }
    data->partial_sums[rank] = sum;
    synchronize();

    for (int step = 1; step < NUM_THREADS; step *= 2) {
        if ((rank % (2 * step)) == (2 * step - 1)) {
            data->partial_sums[rank] += data->partial_sums[rank - step];
        }
        synchronize();
    }
}

typedef struct {
    doublereal *val;
    int *col_ind;
    int *col_ptr;
    doublereal *x;
    doublereal *result;
} matvec_args;

matvec_args *mv_arg;

void matvec(int rank, void* arg) {
    matvec_args *data = (matvec_args*) arg;
    int small = n / NUM_THREADS;
    int large = small + 1;
    int num_large = n % NUM_THREADS;
    int block = (rank < num_large ? large : small);
    int first = (rank < num_large) ? rank * large : num_large * large + (rank - num_large) * small;
    int last = first + block;

    for (int i = first; i < last; i++) {
        data->result[i] = 0.0;
        for (int j = data->col_ptr[i]; j < data->col_ptr[i + 1]; j++) {
            data->result[i] += data->val[j] * data->x[data->col_ind[j]];
        }
    }
}

// --- Conjugate Gradient Method ---
void conjugate_gradient(doublereal *A, doublereal *b, doublereal *x) {
    doublereal *r = calloc(n, sizeof(doublereal));
    doublereal *d = calloc(n, sizeof(doublereal));
    doublereal *Ad = calloc(n, sizeof(doublereal));
    doublereal *tmp = calloc(n, sizeof(doublereal));
    doublereal alpha_k, beta_k, r_dot, r_new_dot, b_norm;

    // Convert A to CSR format
    int nonZero = 0;
    for (int i = 0; i < n * n; i++) if (A[i] != 0) nonZero++;

    mv_arg = malloc(sizeof(matvec_args));
    mv_arg->val = malloc(nonZero * sizeof(doublereal));
    mv_arg->col_ind = malloc(nonZero * sizeof(int));
    mv_arg->col_ptr = calloc(n + 1, sizeof(int));
    mv_arg->x = x;
    mv_arg->result = calloc(n, sizeof(doublereal));

    int idx = 0;
    for (int i = 0; i < n; i++) {
        mv_arg->col_ptr[i] = idx;
        for (int j = 0; j < n; j++) {
            if (A[j * n + i] != 0) {
                mv_arg->val[idx] = A[j * n + i];
                mv_arg->col_ind[idx++] = j;
            }
        }
    }
    mv_arg->col_ptr[n] = idx;

    // r = b - A * x
    function = matvec;
    argument = mv_arg;
    synchronize();
    matvec(0, argument);
    synchronize();

    vec_op_arg = malloc(sizeof(vec_op_args));
    vec_op_arg->a = b;
    vec_op_arg->b = mv_arg->result;
    vec_op_arg->alpha = 1.0;
    vec_op_arg->beta = -1.0;
    vec_op_arg->result = tmp;
    function = vec_op;
    argument = vec_op_arg;
    synchronize();
    vec_op(0, argument);
    synchronize();

    vec_copy_arg = malloc(sizeof(vec_copy_args));
    vec_copy_arg->src = tmp;
    vec_copy_arg->dest = r;
    function = vec_copy;
    argument = vec_copy_arg;
    synchronize();
    vec_copy(0, argument);
    synchronize();

    vec_copy_arg->src = r;
    vec_copy_arg->dest = d;
    synchronize();
    vec_copy(0, argument);
    synchronize();

    dot_arg = malloc(sizeof(dot_args));
    dot_arg->a = r;
    dot_arg->b = r;
    dot_arg->partial_sums = calloc(NUM_THREADS, sizeof(doublereal));
    function = dot_product;
    argument = dot_arg;
    synchronize();
    dot_product(0, argument);
    synchronize();
    r_dot = dot_arg->partial_sums[NUM_THREADS - 1];

    dot_arg->a = b;
    dot_arg->b = b;
    synchronize();
    dot_product(0, argument);
    synchronize();
    b_norm = sqrt(dot_arg->partial_sums[NUM_THREADS - 1]);

    while (sqrt(r_dot) / b_norm > TOLERANCE) {
        // Ad = A * d
        mv_arg->x = d;
        function = matvec;
        argument = mv_arg;
        synchronize();
        matvec(0, argument);
        synchronize();

        vec_copy_arg->src = mv_arg->result;
        vec_copy_arg->dest = Ad;
        function = vec_copy;
        argument = vec_copy_arg;
        synchronize();
        vec_copy(0, argument);
        synchronize();

        dot_arg->a = d;
        dot_arg->b = Ad;
        function = dot_product;
        argument = dot_arg;
        synchronize();
        dot_product(0, argument);
        synchronize();
        alpha_k = r_dot / dot_arg->partial_sums[NUM_THREADS - 1];

        // x = x + alpha * d
        vec_op_arg->a = x;
        vec_op_arg->b = d;
        vec_op_arg->alpha = 1.0;
        vec_op_arg->beta = alpha_k;
        vec_op_arg->result = tmp;
        function = vec_op;
        argument = vec_op_arg;
        synchronize();
        vec_op(0, argument);
        synchronize();

        vec_copy_arg->src = tmp;
        vec_copy_arg->dest = x;
        synchronize();
        vec_copy(0, argument);
        synchronize();

        // r = r - alpha * Ad
        vec_op_arg->a = r;
        vec_op_arg->b = Ad;
        vec_op_arg->alpha = 1.0;
        vec_op_arg->beta = -alpha_k;
        vec_op_arg->result = tmp;
        function = vec_op;
        argument = vec_op_arg;
        synchronize();
        vec_op(0, argument);
        synchronize();

        vec_copy_arg->src = tmp;
        vec_copy_arg->dest = r;
        synchronize();
        vec_copy(0, argument);
        synchronize();

        dot_arg->a = r;
        dot_arg->b = r;
        synchronize();
        dot_product(0, argument);
        synchronize();
        r_new_dot = dot_arg->partial_sums[NUM_THREADS - 1];

        beta_k = r_new_dot / r_dot;

        // d = r + beta * d
        vec_op_arg->a = r;
        vec_op_arg->b = d;
        vec_op_arg->alpha = 1.0;
        vec_op_arg->beta = beta_k;
        vec_op_arg->result = tmp;
        synchronize();
        vec_op(0, argument);
        synchronize();

        vec_copy_arg->src = tmp;
        vec_copy_arg->dest = d;
        synchronize();
        vec_copy(0, argument);
        synchronize();

        r_dot = r_new_dot;
    }

    free(r); free(d); free(Ad); free(tmp);
}

int main() {
    int p = numThreads - 1;
    int max_num_threads = p + 1;
    pthread_t *threads = malloc(p * sizeof(pthread_t));
    wd = malloc(p * sizeof(struct thread_data*));
    for (int i = 0; i < p; i++) wd[i] = malloc(sizeof(struct thread_data));

    pthread_mutex_init(&mutex, 0);
    pthread_cond_init(&cond, 0);
    limit = max_num_threads;
    counter = 0;

    stopped = 0;
    for (int i = 0; i < p; i++) {
        wd[i]->rank = i + 1;
        pthread_create(&threads[i], 0, waiting_thread_function, (void*) wd[i]);
    }

    doublereal *T = malloc(n * n * sizeof(doublereal));
    // Matrix T construction omitted for brevity

    doublereal *y = calloc(n, sizeof(doublereal));
    doublereal *exactSolution = malloc(n * sizeof(doublereal));
    dlaset_("A", &n, &intOne, &one, &one, exactSolution, &n);

    cg(n, T, exactSolution, y, 1e-8);

    integer *ipiv = malloc(n * sizeof(integer));
    dgesv_(&n, &intOne, T, &n, ipiv, exactSolution, &n, &INFO);

    doublereal *residual = malloc(n * sizeof(doublereal));
    for (int i = 0; i < n; i++)
        residual[i] = y[i] - exactSolution[i];

    printf("Relative error (norm): %g\n", dnrm2_(&n, residual, &intOne) / dnrm2_(&n, exactSolution, &intOne));

    // Cleanup and thread termination omitted for brevity

    return 0;
}
