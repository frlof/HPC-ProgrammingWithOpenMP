# High Performance Computing

This repository contains assignment `Programming with OpenMP` and `MPI Programming` from the course *DD2356* at *KTH*.

**Group:** 3

**Group Members:** Fredrik Löf, Victor Carle

## MPI

All the files for this assignment are in the `MPI` folder.

### Build

Instructions for running the code on `Beskow`.

1. *Start by setting up the environment:*
   1. **Kinit** to get ssh key to beskow  
      *Example: `kinit --forward username@NADA.KTH.SE`*
   2. **ssh** into beskow  
      *Example: `ssh username@beskow.pdc.kth.se`*
   3. **salloc** to allocate node  
      *Example: `salloc -t 1:00:00 -A edu19.DD2356 --nodes=2`*
   4. **module  swap** to change compiler to gnu  
      *Example: `module  swap  PrgEnv-cray  PrgEnv-gnu`*
   5. **module load** git  
      *Example: `module load git`*
   6. **git clone** to download the code from this repository  
      *Example: `git clone https://github.com/frlof/HPC-ProgrammingWithOpenMP.git`*

2. *Run the code:*
   1. **make** to comile the code  
      *NOTE: it is important to be in the exercise_x directory when running this command*  
      *Example: `make`*
   2. **aprun** to run the code  
      *NOTE: flag `n` specifies the number of threads, flag `N` specifies the number of cores*  
      *Example: `aprun -n 64 -N 32 DirectoryToBinary`*

**NOTE: Every program needs its specified input to run properly**

### Exercise 6
`verify.c` in the `dat` folder is used to check whether the solution is correct.

#### Results

|      |      |      |      |      |
| ---- | ---- | ---- | ---- | ---- |
|      |      |      |      |      |
|      |      |      |      |      |
|      |      |      |      |      |
