# Parallel Merge Sort (Full Threaded Version)

This project implements a fully threaded **parallel merge sort** using **POSIX threads (pthreads)**. It sorts an array of floating-point numbers by:
1. Sorting subarrays in parallel,
2. Merging them using a tree-structured recursive merge,
3. Synchronizing thread execution using custom barriers.

## Files

- `parallel_merge_sort_full.c` – Main source file
- `HelpFunctions.h` – Must include:
  - `void sequential_merge_sort(double *a, int p, int r, double *b)`
  - `void parallel_copy(int rank, void *arg)`
  - `void merge_recursion(int rank, void *arg)`

## Highlights

- Custom thread barrier using mutex and condition variables
- General-purpose thread pool executing assigned functions
- Multiple levels of parallel merging with thread reuse
- Fully dynamic allocation of thread work ranges

## Compilation

```bash
gcc -o parallel_sort parallel_merge_sort_full.c -pthread
