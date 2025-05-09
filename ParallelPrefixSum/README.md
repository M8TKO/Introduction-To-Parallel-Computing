# Parallel Prefix Sum (Scan) using Pthreads

This project implements a **parallel prefix sum algorithm** (also known as a scan operation) using **POSIX threads (pthreads)** in C.

Each thread computes the prefix sum for a chunk of the input array, and the results are synchronized using barriers to ensure correctness during the upward and downward sweep phases.

## Files

- `parallel_prefix_sum.c` – Main source file containing the implementation.

## How It Works

The prefix sum operation is performed in multiple phases:
1. **Local scan:** Each thread calculates the partial sum of its assigned block.
2. **Upward sweep:** Partial results are merged in a tree-like fashion.
3. **Downward sweep:** Final prefix sums are propagated back to all blocks.
4. **Final assembly:** Each thread fills in its segment of the result array.

Synchronization is handled using a custom barrier implemented with `pthread_mutex_t` and `pthread_cond_t`.

## Compilation

Use `gcc` and link with the `-pthread` flag:

```bash
gcc -o prefix_sum parallel_prefix_sum.c -pthread
