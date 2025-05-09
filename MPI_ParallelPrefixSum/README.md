# MPI Parallel Prefix Sum (Scan)

This project implements a **parallel prefix sum (scan)** using **MPI (Message Passing Interface)**. The computation is structured as a tree where:

- Leaf processes compute local prefix sums of blocks,
- Internal nodes aggregate and propagate sums up and down the tree.

## Files

- `mpi_parallel_prefix_sum.c` – Main implementation using MPI.

## Algorithm Overview

The parallel prefix sum is split into two major phases:
1. **Upward sweep (reduction):** Each leaf computes a local sum and sends it up the tree.
2. **Downward sweep (distribution):** Root distributes prefix results back to leaves to complete the scan.

## Build Instructions

Make sure you have `mpicc` installed (from OpenMPI or MPICH):

```bash
mpicc -o mpi_scan mpi_parallel_prefix_sum.c
