#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include "HelpFunctions.h"

// Synchronization structure
typedef struct synchronize_structure {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int counter;
    int limit;
} synchronize_type;

synchronize_type *synchronize_all_vars;

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

// Function pointer for thread work
typedef void (*working_function_t)(int, void*);
working_function_t function;

static void* argument;

// Thread structure
struct thread_data {
    int rank;
};

struct thread_data **wd;
int stopped;

// Worker thread loop
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

// Parallel sort
typedef struct {
    int num_threads;
    double *a;
    int n;
    int *borders;
    double *b;
} parallel_sort_type;

parallel_sort_type *parallel_sort_arg;

void parallel_sort(int rank, void *arg) {
    parallel_sort_type *pca = (parallel_sort_type*) arg;
    int small = pca->n / pca->num_threads;
    int large = small + 1;
    int num_large = pca->n % pca->num_threads;
    int first, last;
    int myp = rank;
    int block = (myp < num_large ? large : small);

    if (myp < num_large)
        first = myp * large;
    else
        first = num_large * large + (myp - num_large) * small;
    last = first + block;

    pca->borders[rank] = first;
    if (last == pca->n)
        pca->borders[pca->num_threads] = pca->n;

    sequential_merge_sort(pca->a, first, last - 1, pca->b);
}

// Parallel merge
typedef struct {
    int num_threads;
    double *a;
    int *borders;
    double *b;
    synchronize_type **synchronize_nt_vars;
} multiple_parallel_merge_type;

multiple_parallel_merge_type *multiple_parallel_merge_arg;

merge_recursion_type **merge_recursion_arg;
parallel_copy_type **parallel_copy_arg;

void multiple_parallel_merge(int rank, void *arg) {
    multiple_parallel_merge_type *mpma = (multiple_parallel_merge_type*) arg;
    int merge_num = rank / mpma->num_threads;
    int step = mpma->num_threads / 2;

    merge_recursion_arg[merge_num]->a = mpma->a;
    merge_recursion_arg[merge_num]->b = mpma->b;
    merge_recursion_arg[merge_num]->num_threads = mpma->num_threads;
    merge_recursion_arg[merge_num]->p1 = mpma->borders[2 * merge_num * step];
    merge_recursion_arg[merge_num]->r1 = mpma->borders[2 * merge_num * step + step] - 1;
    merge_recursion_arg[merge_num]->p2 = mpma->borders[2 * merge_num * step + step];
    merge_recursion_arg[merge_num]->r2 = mpma->borders[2 * merge_num * step + 2 * step] - 1;
    merge_recursion_arg[merge_num]->p3 = merge_recursion_arg[merge_num]->p1;

    merge_recursion(rank % mpma->num_threads, merge_recursion_arg[merge_num]);
    synchronize(mpma->synchronize_nt_vars[merge_num]);

    parallel_copy_arg[merge_num]->a = mpma->a + merge_recursion_arg[merge_num]->p1;
    parallel_copy_arg[merge_num]->b = mpma->b + merge_recursion_arg[merge_num]->p1;
    parallel_copy_arg[merge_num]->num_threads = mpma->num_threads;
    parallel_copy_arg[merge_num]->n = merge_recursion_arg[merge_num]->r2 - merge_recursion_arg[merge_num]->p1 + 1;

    parallel_copy(rank % mpma->num_threads, parallel_copy_arg[merge_num]);
}

void parallel_merge_sort(double *a, int n) {
    int max_num_threads = parallel_sort_arg->num_threads;
    int h = log2(max_num_threads);
    int m = max_num_threads / 2;
    int nt_merge = 2;

    int *borders = malloc((max_num_threads + 1) * sizeof(int));
    double *b = malloc(n * sizeof(double));

    parallel_sort_arg->b = b;
    parallel_sort_arg->borders = borders;

    function = parallel_sort;
    argument = parallel_sort_arg;
    synchronize(synchronize_all_vars);
    parallel_sort(0, argument);
    synchronize(synchronize_all_vars);

    multiple_parallel_merge_arg->borders = borders;
    multiple_parallel_merge_arg->b = b;

    for (int k = 1; k <= h; k++) {
        for (int i = 0; i < m; i++)
            multiple_parallel_merge_arg->synchronize_nt_vars[i]->limit = nt_merge;

        multiple_parallel_merge_arg->num_threads = nt_merge;
        function = multiple_parallel_merge;
        argument = multiple_parallel_merge_arg;
        synchronize(synchronize_all_vars);
        multiple_parallel_merge(0, argument);
        synchronize(synchronize_all_vars);

        m /= 2;
        nt_merge *= 2;
    }

    free(b);
    free(borders);
}

int main() {
    int n = 100;
    double a[] = {99, 32, 40, 22, 34, 92, 91, 35, 6, 55,
                  3, 96, 68, 16, 69, 11, 54, 30, 45, 77,
                  60, 74, 78, 72, 62, 70, 51, 33, 7, 86,
                  38, 100, 58, 76, 81, 89, 42, 28, 17, 41,
                  47, 98, 80, 14, 46, 56, 63, 93, 8, 67,
                  84, 90, 97, 83, 59, 79, 5, 48, 53, 29,
                  21, 25, 52, 37, 64, 31, 49, 27, 61, 88,
                  50, 87, 26, 43, 94, 19, 44, 15, 73, 1,
                  36, 82, 71, 23, 65, 2, 4, 18, 85, 75,
                  24, 95, 39, 13, 9, 66, 20, 57, 10, 12};

    int max_num_threads = 16;

    synchronize_all_vars = malloc(sizeof(synchronize_type));
    pthread_mutex_init(&synchronize_all_vars->mutex, 0);
    pthread_cond_init(&synchronize_all_vars->cond, 0);
    synchronize_all_vars->limit = max_num_threads;
    synchronize_all_vars->counter = 0;

    pthread_t *threads = malloc((max_num_threads - 1) * sizeof(pthread_t));
    wd = malloc((max_num_threads - 1) * sizeof(struct thread_data*));
    for (int i = 0; i < max_num_threads - 1; i++)
        wd[i] = malloc(sizeof(struct thread_data));

    stopped = 0;
    for (int i = 0; i < max_num_threads - 1; i++) {
        wd[i]->rank = i + 1;
        pthread_create(&threads[i], 0, waiting_thread_function, (void*) wd[i]);
    }

    int pom = max_num_threads / 2;
    merge_recursion_arg = malloc(pom * sizeof(merge_recursion_type*));
    parallel_copy_arg = malloc(pom * sizeof(parallel_copy_type*));

    for (int i = 0; i < pom; i++) {
        merge_recursion_arg[i] = malloc(sizeof(merge_recursion_type));
        parallel_copy_arg[i] = malloc(sizeof(parallel_copy_type));
    }

    parallel_sort_arg = malloc(sizeof(parallel_sort_type));
    parallel_sort_arg->num_threads = max_num_threads;
    parallel_sort_arg->a = a;
    parallel_sort_arg->n = n;

    multiple_parallel_merge_arg = malloc(sizeof(multiple_parallel_merge_type));
    multiple_parallel_merge_arg->synchronize_nt_vars = malloc(pom * sizeof(synchronize_type*));
    multiple_parallel_merge_arg->a = a;
    multiple_parallel_merge_arg->num_threads = max_num_threads;

    for (int i = 0; i < pom; i++) {
        multiple_parallel_merge_arg->synchronize_nt_vars[i] = malloc(sizeof(synchronize_type));
        multiple_parallel_merge_arg->synchronize_nt_vars[i]->counter = 0;
        pthread_mutex_init(&multiple_parallel_merge_arg->synchronize_nt_vars[i]->mutex, 0);
        pthread_cond_init(&multiple_parallel_merge_arg->synchronize_nt_vars[i]->cond, 0);
    }

    parallel_merge_sort(a, n);

    for (int i = 0; i < n; i++)
        printf("%g\n", a[i]);

    stopped = 1;
    synchronize(synchronize_all_vars);

    for (int i = 0; i < max_num_threads - 1; i++)
        pthread_join(threads[i], 0);

    for (int i = 0; i < max_num_threads - 1; i++)
        free(wd[i]);

    free(wd);
    free(threads);
    pthread_mutex_destroy(&synchronize_all_vars->mutex);
    pthread_cond_destroy(&synchronize_all_vars->cond);
    free(synchronize_all_vars);

    for (int i = 0; i < pom; i++) {
        free(merge_recursion_arg[i]);
        free(parallel_copy_arg[i]);
    }
    free(merge_recursion_arg);
    free(parallel_copy_arg);
    free(parallel_sort_arg);

    for (int i = 0; i < pom; i++) {
        pthread_mutex_destroy(&multiple_parallel_merge_arg->synchronize_nt_vars[i]->mutex);
        pthread_cond_destroy(&multiple_parallel_merge_arg->synchronize_nt_vars[i]->cond);
        free(multiple_parallel_merge_arg->synchronize_nt_vars[i]);
    }
    free(multiple_parallel_merge_arg->synchronize_nt_vars);
    free(multiple_parallel_merge_arg);

    return 0;
}
