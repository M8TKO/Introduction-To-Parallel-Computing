#ifndef HELPFUNCTIONS
#define HELPFUNCTIONS

int binary_search(double *a, int p, int r, double x);
void sequential_merge(double *a, int p1,int r1, int p2, int r2, double *b, int p3);

typedef struct merge_recursion_structure{
    int num_threads;
    double *a;
    int p1;
    int r1;
    int p2;
    int r2;
    double *b;
    int p3;
} merge_recursion_type;
merge_recursion_type **merge_recursion_arg;
void merge_recursion( int rank,merge_recursion_type *arg);

typedef struct parallel_copy_structure{
    int num_threads;
    int n;
    double *a;
    double *b;
} parallel_copy_type;
parallel_copy_type **parallel_copy_arg;
void parallel_copy( int rank, parallel_copy_type *arg);

void sequential_merge_sort(double *a, int p, int r, double *b);

#endif 