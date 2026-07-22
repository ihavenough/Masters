# C Finite Element Method (FEM) Solver

**Authors:** Karolina Ochman, Beining San  
**Date:** May 2026  

## Overview
This repository contains a high-performance Finite Element Method (FEM) program written in C (with Fortran solvers). The solver is designed to numerically resolve elliptic partial differential equations (PDEs) using **P1 elements on triangles** and **Q1 elements on quadrangles**. 

The program was tested on multiple domains (a 1x1 square and a 1/3x1/3 square) and handles various boundary conditions, including:
* Homogeneous Dirichlet
* Non-homogeneous Dirichlet
* Neumann

The primary objective of this project is to study and validate the order of convergence of the method.

## Project Structure
The codebase strictly separates the core mathematical logic from the physical definitions:

* `main.c` / `mainerr.c`: The core entry points. They read the mesh, trigger the assembly, manage matrix conversion, execute the Cholesky resolution, and calculate the final errors.
* `tp2b.c` (Physics Dictionary): Determines the physical configuration based on the selected test case (`nucas`). It provides the exact solution (`solex`), the right-hand side source term (`fomega`), and defines whether a boundary is Dirichlet or Neumann.
* `cal1Elem.c` (Local Mathematical Engine): Computes surface and edge integrals on the reference elements and enforces Dirichlet boundary conditions.
* `assemblage.c` & `dsmoapr.c`: Manages the construction of the global matrix (S.M.D, S.M.O) and its conversion into a Skyline/Profile envelope structure required by the Fortran solver.

## How to Build and Run

### 1. Interactive Mode
In this mode, the program prompts the user to select a test case (1-6) and input the name of a specific mesh file (e.g., `d1t1_2`).

**Compilation:**
```bash
gcc main.c assemblage.c lecture.c utils.c cal1Elem.c tp2b.c dsmoapr.c assmat.f affsmd.f cdesse.f affsmo.f tri.f ltlpr.f rsprl.f rspru.f affsol.f impmpr.f -o prog -lm -lgfortran

**Execution:**
```bash
./prog

### 2. Batch Mode (Automated Convergence Analysis)
To automate the convergence analysis, we developed a dedicated batch pipeline (mainerr.c). This script loops through all mathematical test cases (1 to 6) across their respective mesh sequences, from the coarsest to the finest.
It automatically assembles and solves the linear system for each configuration, computes the relative quadratic errors (L2) and relative maximum errors (L∞), and outputs the progress to the console.

**Compilation:**
```bash
gcc mainerr.c assemblage.c lecture.c utils.c cal1Elem.c tp2b.c tp4.c dsmoapr.c solex.c assmat.f affsmd.f cdesse.f affsmo.f tri.f ltlpr.f rsprl.f rspru.f affsol.f impmpr.f -o prog -lm -lgfortran

**Execution:**
```bash
./prog

Note: This mode requires no interactive input. Upon completion, all error data is automatically extracted and saved to erreurs.txt, which can then be imported into Python to plot the convergence curves.
