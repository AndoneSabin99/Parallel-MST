# Parallel-MST

University project for HPC4DS course

This project implements the Kruskal's algorithm for creating a minimum spanning tree (MST) of a weighted, undirected graph in C language with parallelization via mpi and openmp.

Kruskal's algorithm works as following: first sort all the edges via parallelized sorting algorithm (for this project, merge sort has been choosen as the sorting algorithm); then, keep track of the components via union-find data structures with union by rank and path compression.

## Commands to compile

First compile the .c files like following:

gcc -c -g data_structures.c -o data_structures.o

mpicc -std=gnu99 -g -Wall -fopenmp -c sort.c -o sort.o

mpicc -std=gnu99 -g -Wall -fopenmp -c mst.c -o mst.o



Then create the library files for the sort and data_structure files:

ar rcs libsort.a sort.o

ar rcs libdata_structures.a data_structures.o



Then, once we have compiled all the files, create the executable file:

mpicc -std=gnu99 -g -Wall -fopenmp -o mst mst.o data_structures.o sort.o -L<path_libraries> -lsort -ldata_structures


## Commands to run (inside an hpc cluster that works with jobs):

qsub run.sh

To run directly (not recommended to do this directly on the login node of the hpc cluster, but rather write this inside the .sh file)

mpirun.actual -np 4 <complete_path/mst> <complete_path/datasets/input_file.txt> 4
