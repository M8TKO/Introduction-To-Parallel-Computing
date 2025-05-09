#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include "HelpFunctions.h"

// Synchronization structure and global variables
typedef struct synchronize_structure {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int counter;
    int limit;
} synchronize_type;

synchronize_type *synchronize_all_vars, **synchronize_nt_vars;

void synchronize(synchronize_type *synchronize_vars) {
    pthread_mutex_lock(&synchronize_vars->mutex);
    synchronize_vars->counter++;

    if (synchronize_vars->counter == synchronize_vars->limit) {
        synchronize_vars->counter = 0;
        pthread_cond_broadcast(&synchronize_vars->cond);
    } else {
        pthread_cond_wait(&synchronize_vars->cond, &synchronize_vars->mutex);
    }

    pthread_mutex_unlock(&synchronize_vars->mutex);
}

// Function pointer for work assignment
typedef void (*working_function_t)(int, void*);
working_function_t function;

static void* argument;

// Thread data structure
struct thread_data {
    int rank;
};

struct thread_data **wd;
int stopped;

// Thread loop function
static void *waiting_thread_function(void *arg) {
    struct thread_data *thread_data = (struct thread_data*) arg;
    const int rank = thread_data->rank;

    for (;;) {
        synchronize(synchronize_all_vars);
        if (stopped)
            break;
        else {
            function(rank, argument);
            synchronize(synchronize_all_vars);
        }
    }
    return NULL;
}

// Argument for merge recursion
typedef struct merge_recursion_structure {
    int num_threads;
    double *a;
    int p1, r1;
    int p2, r2;
    double *b;
    int p3;
} merge_recursion_type;

merge_recursion_type *merge_recursion_arg;

void merge_recursion(int rank, void *arg) {
    merge_recursion_type *mra = (merge_recursion_type*) arg;
    int nt = mra->num_threads, p1 = mra->p1, p2 = mra->p2,
        r1 = mra->r1, r2 = mra->r2, p3 = mra->p3;
    int k = log2(nt);
    int break_point = nt / 2;
    int step = break_point / 2;

    int q1, q2, q3;
    double x;
    for (int i = 0; i < k; i++) {
        // Split points q1 & q2
        q1 = (p1 + r1) / 2;
        x = mra->a[q1];
        q2 = binary_search(mra->a, p2, r2, x);
        q3 = p3 + (q1 - p1) + (q2 - p2);
        mra->b[q3] = x;

        if (rank < break_point) {
            break_point -= step;
            r1 = q1 - 1;
            r2 = q2 - 1;
        } else {
            break_point += step;
            p1 = q1 + 1;
            p2 = q2;
            p3 = q3 + 1;
        }
        step /= 2;
    }

    sequential_merge(mra->a, p1, r1, p2, r2, mra->b, p3);
}

// Argument for parallel copy
typedef struct parallel_copy_structure {
    int num_threads;
    int n;
    double *a;
    double *b;
} parallel_copy_type;

parallel_copy_type *parallel_copy_arg;

void parallel_copy(int rank, void *arg) {
    parallel_copy_type *pca = (parallel_copy_type*) arg;
    int small = pca->n / pca->num_threads;
    int large = small + 1;
    int num_large = pca->n % pca->num_threads;
    int first, last;
    int myp = rank;
    int block = (myp < num_large ? large : small);

    if (myp < num_large) {
        first = myp * large;
    } else {
        first = num_large * large + (myp - num_large) * small;
    }
    last = first + block;

    for (int i = 0; i < last - first; i++) {
        pca->a[first + i] = pca->b[first + i];
    }
}

void parallel_merge(double *a, int p, int q, int r) {
    double* b = malloc((r - p + 1) * sizeof(double));

    // Merge recursion
    merge_recursion_arg->b = b;
    function = merge_recursion;
    argument = merge_recursion_arg;

    synchronize(synchronize_all_vars);
    merge_recursion(0, argument);
    synchronize(synchronize_all_vars);

    // Parallel copy
    parallel_copy_arg->b = b;
    function = parallel_copy;
    argument = parallel_copy_arg;

    synchronize(synchronize_all_vars);
    parallel_copy(0, argument);
    synchronize(synchronize_all_vars);

    free(b);
}

int main() {
    int n = 100;
    double a[] = {
        2, 4, 5, 6, 7, 11, 14, 16, 17, 18,
        19, 20, 21, 22, 23, 26, 27, 28, 30, 31,
        33, 34, 35, 36, 37, 41, 42, 43, 45, 46,
        47, 51, 53, 55, 57, 58, 61, 64, 65, 66,
        67, 68, 69, 73, 76, 82, 83, 84, 85, 87,
        88, 89, 90, 91, 92, 94, 95, 98, 99, 100,
        1, 3, 8, 9, 10, 12, 13, 15, 24, 25,
        29, 32, 38, 39, 40, 44, 48, 49, 50, 52,
        54, 56, 59, 60, 62, 63, 70, 71, 72, 74,
        75, 77, 78, 79, 80, 81, 86, 93, 96, 97
    };

    int p = 15;
    int max_num_threads = p + 1;

    synchronize_all_vars = malloc(sizeof(synchronize_type));
    pthread_mutex_init(&synchronize_all_vars->mutex, 0);
    pthread_cond_init(&synchronize_all_vars->cond, 0);
    synchronize_all_vars->limit = max_num_threads;
    synchronize_all_vars->counter = 0;

    pthread_t *threads = malloc(p * sizeof(pthread_t));
    wd = malloc(p * sizeof(struct thread_data*));

    for (int i = 0; i < p; i++)
        wd[i] = malloc(sizeof(struct thread_data));

    // Start worker threads
    stopped = 0;
    for (int i = 0; i < p; i++) {
        wd[i]->rank = i + 1;
        pthread_create(&threads[i], 0, waiting_thread_function, (void*) wd[i]);
    }

    // Set up merge recursion arguments
    merge_recursion_arg = malloc(sizeof(merge_recursion_type));
    merge_recursion_arg->num_threads = max_num_threads;
    merge_recursion_arg->a = a;
    merge_recursion_arg->p1 = 0;
    merge_recursion_arg->r1 = 59;
    merge_recursion_arg->p2 = 60;
    merge_recursion_arg->r2 = 99;
    merge_recursion_arg->p3 = 0;

    // Set up parallel copy arguments
    parallel_copy_arg = malloc(sizeof(parallel_copy_type));
    parallel_copy_arg->n = n;
    parallel_copy_arg->num_threads = max_num_threads;
    parallel_copy_arg->a = a;

    // Perform parallel merge and print result
    parallel_merge(a, 0, 59, 99);

    for (int i = 0; i < n; i++)
        printf("%g\n", a[i]);

    // Stop threads
    stopped = 1;
    synchronize(synchronize_all_vars);

    for (int i = 0; i < p; i++)
        pthread_join(threads[i], 0);

    for (int i = 0; i < p; i++)
        free(wd[i]);

    free(wd);
    free(threads);
    free(merge_recursion_arg);
    free(parallel_copy_arg);
    pthread_mutex_destroy(&synchronize_all_vars->mutex);
    pthread_cond_destroy(&synchronize_all_vars->cond);
    free(synchronize_all_vars);

    return 0;
}
