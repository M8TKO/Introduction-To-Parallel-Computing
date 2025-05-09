# Parallel Merge with Pthreads

This project implements a parallel merge operation using **POSIX threads (pthreads)**. It merges two sorted halves of an array into a final sorted array using multiple threads working in coordination.

## Files

- `parallel_merge_sort_threads.c` – Main implementation using pthreads
- `HelpFunctions.h` – Header file containing `sequential_merge()` and `binary_search()` (assumed provided)

## Algorithm Overview

The implementation consists of two main steps:

1. **Recursive Parallel Merge**
   - Threads are used to simulate levels of a recursive binary merge.
   - The merge splits the array, searches for element positions, and builds a merged result in a buffer.

2. **Parallel Copy**
   - The merged buffer is copied back into the original array in parallel.

Threads are synchronized via a custom `synchronize()` function that mimics a barrier.

## Requirements

- POSIX threads (`pthread.h`)
- A header file `HelpFunctions.h` defining:
  - `int binary_search(double *array, int left, int right, double key);`
  - `void sequential_merge(double *a, int p1, int r1, int p2, int r2, double *b, int p3);`

## Compilation

```bash
gcc -o parallel_merge parallel_merge_sort_threads.c -pthread
