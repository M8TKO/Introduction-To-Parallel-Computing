# Parallel Computing – Course Repository

This repository contains all five programming assignments and the final project completed by **Matko Petričić** for the **Parallel Computing** course.

Each assignment focuses on a fundamental parallel programming concept and is implemented in **C**, using **Pthreads**, **MPI**, or **BLAS/LAPACK**. The final project is a complete parallel solver for linear systems using the **Conjugate Gradient method**.

---

## 📘 What I Learned

Throughout this course, I gained hands-on experience in:

- 🧵 Writing and synchronizing multithreaded programs with **POSIX threads**
- 📡 Building distributed memory programs using **MPI (Message Passing Interface)**
- ⚙️ Implementing **custom barriers** and thread synchronization
- 🧮 Efficient use of **BLAS** and **LAPACK** for matrix and vector operations
- 📚 Sparse matrix storage (CSR format) and **matrix-vector multiplication**
- 🔄 Recursive and tree-based algorithms for scan and merge
- 📈 Benchmarking and validating results with **serial vs. parallel implementations**
- 🧪 Using **norms** to measure numerical accuracy of parallel algorithms

---

## 🧠 Projects Overview

### Task 1: Parallel Prefix Sum (Pthreads)
Implemented a parallel prefix sum (scan) using Pthreads. Threads compute local prefix sums and propagate results via a tree-based reduction.  
**Keywords:** Barrier sync, thread-local computation

### Task 2: Parallel Prefix Sum (MPI)
Recreated the same scan operation using MPI and a logical binary tree of processes.  
**Keywords:** Scatter/Gather, collective communication

### Task 3: Parallel Merge of Two Sorted Segments
Each thread finds the correct position of its element via binary search and fills a result buffer in parallel.  
**Keywords:** Binary search, divide and conquer

### Task 4: Full Parallel Merge Sort
Merged task 3 recursively into a complete merge sort with multistage merging and re-used threads across stages.  
**Keywords:** Multi-way merging, thread reuse

### Task 5: Blocked Matrix Multiplication (MPI)
Performed parallel matrix multiplication by dividing matrices into column blocks. Verified correctness against LAPACK’s `dgemm`.  
**Keywords:** Block distribution, MPI_Scatterv/Gatherv

---

## 🔴 Final Project: Parallel Conjugate Gradient Solver

Built a complete solver for `Ax = b` where `A` is a sparse SPD matrix:

- Matrix stored in **CSR format**
- Matrix-vector products computed in parallel
- Scalar products and vector updates are distributed across threads
- Validated against LAPACK's direct solver (`dgesv`) using norm-based residual comparison

This project combined everything from the course — synchronization, threading, vector ops, and numerical accuracy — into one cohesive application.

---

## 🛠️ Build Instructions

Each file can be compiled using `gcc` or `mpicc`, depending on whether it's a threaded or MPI program.

Examples:
```bash
# Pthreads-based
gcc -o task1 task1_parallelPrefixSum.c -pthread

# MPI-based
mpicc -o task5 task5_mpiBlockedMatrixMultiplication.c -llapack -lblas

# Final project
gcc -o cgSolver final_parallelConjugateGradientSolver.c -pthread -llapack -lblas
