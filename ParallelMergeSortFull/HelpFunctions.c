#include "HelpFunctions.h"
#include <stdio.h>
#include <math.h>

// Binary search: returns index of the first element >= x in sorted array a[p..r]
int binary_search(double *a, int p, int r, double x) {
    int low = p, high = r + 1, mid, q = p;
    while (low < high) {
        mid = (low + high) / 2;
        if (x <= a[mid])
            high = mid;
        else
            low = mid + 1;
        q = low;
    }
    return q;
}

// Merges two sorted segments of array a into b starting at index p3
void sequential_merge(double *a, int p1, int r1, int p2, int r2, double *b, int p3) {
    int n1 = r1 - p1 + 1;
    int n2 = r2 - p2 + 1;
    int i = 0, j = 0, k = p3;

    while (i < n1 && j < n2) {
        if (a[p1 + i] <= a[p2 + j])
            b[k++] = a[p1 + i++];
        else
            b[k++] = a[p2 + j++];
    }

    while (i < n1)
        b[k++] = a[p1 + i++];
    while (j < n2)
        b[k++] = a[p2 + j++];
}

// Recursive merge for parallel usage
void merge_recursion(int rank, merge_recursion_type *arg) {
    int nt = arg->num_threads;
    int p1 = arg->p1, r1 = arg->r1, p2 = arg->p2, r2 = arg->r2;
    int p3 = arg->p3;

    if (p1 > r1 && p2 > r2)
        return;
    if ((r1 - p1) < (r2 - p2)) {
        int temp;
        temp = p1; p1 = p2; p2 = temp;
        temp = r1; r1 = r2; r2 = temp;
    }

    int k = log2(nt);
    int break_point = nt / 2;
    int step = break_point / 2;

    int q1, q2, q3;
    double x;
    for (int i = 0; i < k; i++) {
        q1 = (p1 + r1) / 2;
        x = arg->a[q1];
        q2 = binary_search(arg->a, p2, r2, x);
        q3 = p3 + (q1 - p1) + (q2 - p2);
        arg->b[q3] = x;

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

    sequential_merge(arg->a, p1, r1, p2, r2, arg->b, p3);
}

// Copies array b to a in parallel chunks based on thread rank
void parallel_copy(int rank, parallel_copy_type *arg) {
    int n = arg->n;
    int num_threads = arg->num_threads;

    int small = n / num_threads;
    int large = small + 1;
    int num_large = n % num_threads;
    int first, last;
    int myp = rank;
    int block = (myp < num_large ? large : small);

    if (myp < num_large)
        first = myp * large;
    else
        first = num_large * large + (myp - num_large) * small;
    last = first + block;

    for (int i = 0; i < last - first; i++)
        arg->a[first + i] = arg->b[first + i];
}

// Standard sequential merge sort using temporary buffer b
void sequential_merge_sort(double *a, int p, int r, double *b) {
    if (p >= r)
        return;

    int q = (p + r) / 2;
    sequential_merge_sort(a, p, q, b);
    sequential_merge_sort(a, q + 1, r, b);
    sequential_merge(a, p, q, q + 1, r, b, p);

    for (int i = p; i <= r; i++)
        a[i] = b[i];
}
