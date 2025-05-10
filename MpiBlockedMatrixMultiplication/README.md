# MPI Blocked Matrix Multiplication

This project implements **parallel matrix multiplication** using **MPI** with a **blocked and scattered approach**. It splits matrices among processes and performs distributed matrix multiplication using optimized BLAS and LAPACK routines.

## Key Concepts

- **MPI Scatterv/Gatherv**: Flexible distribution of matrix blocks
- **BLAS (dgemm)**: High-performance matrix multiplication
- **CLAPACK**: For matrix initialization and copying
- **Validation**: Result checked against serial `dgemm` multiplication on rank 0

## File

- `mpi_blocked_matrix_multiplication.c` — Main implementation file

## Compilation

Make sure you have MPI, BLAS, and LAPACK installed (e.g., via OpenMPI + Netlib).

```bash
mpicc -o mpi_matmul mpi_blocked_matrix_multiplication.c -llapack -lblas
