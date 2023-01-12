# Parallel-MST

University project for HPC4DS course

This project implements the algorithm of Kruskal for creating a minimum spanning tree (MST) of a weighted, undirected graph in C with parallelization via mpi and openmp.

Kruskal's algorithm works as following: first sort all the edges via parallelized sorting algorithm (for this project, merge sort has been choosen as the sorting algorithm); then, keep track of the components via union-find data structures with union by rank and path compression.

## Commands to Run

To compile:

mpicc -std=gnu99 -g -Wall -fopenmp -o mst mst.c

To run (inside an hpc cluster that works with jobs):

qsub run.sh

To run directly (not recommended to do this directly on the login node of the hpc cluster

mpirun.actual -np 4 <complete_path/mst> <complete_path_input_file> 8
