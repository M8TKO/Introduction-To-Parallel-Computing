# Parallel Conjugate Gradient Solver (C + LAPACK + Pthreads)

This project implements a **parallelized Conjugate Gradient (CG) algorithm** using **POSIX threads** (pthreads) and **LAPACK/BLAS** routines. It solves a linear system `Ax = b` where `A` is a symmetric positive definite matrix, and `b` is known.

## Features

- Parallel matrix-vector multiplication using compressed sparse row (CSR) format
- Parallel scalar products and vector updates using synchronized threads
- Uses LAPACK routines (`dgesv`, `dnrm2`, etc.) for comparison and verification
- Custom thread barrier for synchronization

## File

- `parallelConjugateGradientSolver.c` — Main implementation

## Dependencies

- POSIX threads
- LAPACK and BLAS (e.g., via Netlib or OpenBLAS)
- C compiler (e.g., `gcc`)

## Build Instructions

```bash
gcc -o cg_solver parallelConjugateGradientSolver.c -lpthread -llapack -lblas
