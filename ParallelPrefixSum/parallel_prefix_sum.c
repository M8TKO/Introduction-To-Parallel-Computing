#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
int counter = 0;
int limit;

// Synchronization barrier for threads
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

typedef struct {
    int num_threads;
    int MYPROC;
    int n;
    double *x;
    double *y;
    double *local_y;
} sum_type;

void *sum(void *ptr) {
    sum_type *data = (sum_type*) ptr;
    int nt = data->num_threads;
    int myp = data->MYPROC;
    int dim = data->n;
    double *x = data->x;
    double *y = data->y;
    double *local_y = data->local_y;

    double *local_y2;

    int m = nt;
    int k = log2(nt) + 1;

    // Total number of steps: after the 1st step (i.e., balanced distribution),
    // we have as many numbers as threads. From lectures, we need log2(nt) steps,
    // plus one for the first step.
    int step = 1;

    // Determine the block size and boundaries
    int small = dim / nt;
    int large = small + 1;
    int num_large = dim % nt;
    int first, last;
    int block = (myp < num_large ? large : small);

    if (myp < num_large) {
        first = myp * large;
    } else {
        first = num_large * large + (myp - num_large) * small;
    }
    last = first + block;

    // First step: upward sweep (local prefix sum)
    local_y2 = malloc((last - first) * sizeof(double));
    local_y2[0] = x[first];
    for (int i = first + 1; i < last; i++) {
        local_y2[i - first] = local_y2[i - first - 1] + x[i];
    }

    local_y[myp] = local_y2[last - first - 1];
    synchronize();

    // Continuation of the upward sweep
    m /= 2;
    for (int j = 2; j <= k; j++) {
        for (int i = 1; i <= m; i++) {
            if (myp == i * step - 1) {
                local_y[2 * myp + 1] = local_y[2 * myp - step + 1] + local_y[2 * myp + 1];
            }
        }
        m /= 2;
        step *= 2;
        synchronize();
    }

    // Downward sweep
    m = 1;
    step /= 2;
    for (int j = 1; j <= k - 2; j++) {
        m *= 2;
        step /= 2;

        for (int i = 2; i <= m; i++) {
            if (myp == i * step - 1) {
                local_y[2 * myp - step + 1] = local_y[2 * myp - 2 * step + 1] + local_y[2 * myp - step + 1];
            }
        }
        synchronize();
    }

    // Fill the final y vector
    for (int i = 0; i < nt; i++) {
        if (i == 0) {
            for (int j = first; j < last; j++) {
                y[j] = local_y2[j];
            }
        } else if (myp == i) {
            for (int j = 0; j < last - first; j++) {
                y[first + j] = local_y2[j] + local_y[myp - 1];
            }
        }
    }

    synchronize();
    free(local_y2);
    return NULL;
}

int main() {
    int n = 100;
    int num_threads = 16;
    limit = num_threads;

    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));
    double *local_y = malloc(num_threads * sizeof(double));

    for (int i = 0; i < n; i++) {
        x[i] = 1;
    }

    pthread_t *threads = malloc((num_threads - 1) * sizeof(pthread_t));
    sum_type *arg = malloc(num_threads * sizeof(sum_type));

    for (int i = 0; i < num_threads - 1; i++) {
        arg[i + 1] = (sum_type){num_threads, i + 1, n, x, y, local_y};
        pthread_create(&threads[i], NULL, sum, (void*) (arg + i + 1));
    }

    arg[0] = (sum_type){num_threads, 0, n, x, y, local_y};
    sum((void*) arg);

    for (int i = 0; i < num_threads - 1; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("\nVector local_y:\n");
    for (int i = 0; i < num_threads; i++) {
        printf("%g\n", local_y[i]);
    }

    printf("\nVector y:\n");
    for (int i = 0; i < n; i++) {
        printf("%g\n", y[i]);
    }

    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&cond);

    free(x);
    free(y);
    free(local_y);
    free(arg);
    free(threads);

    return 0;
}
