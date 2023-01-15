# Parallel-MST

University project for HPC4DS course

This project implements the Kruskal's algorithm for creating a minimum spanning tree (MST) of a weighted, undirected graph in C language with parallelization via mpi and openmp.

Kruskal's algorithm works as following: first sort all the edges via parallelized sorting algorithm (for this project, merge sort has been choosen as the sorting algorithm); then, keep track of the components via union-find data structures with union by rank and path compression.

## This reporitory

In this repository it is possible to navigate through two directories:
1) one for the source code <br />
2) one for the dataset, which contains some .txt files with input graphs <br />

Inside the src folder, it is possible to find the source files of the program, together with the libraries and also one example of output of an instance of execution.


## Commands to compile

First compile the .c files like following:

gcc -c -g data_structures.c -o data_structures.o <br />
mpicc -std=gnu99 -g -Wall -fopenmp -c sort.c -o sort.o <br />
mpicc -std=gnu99 -g -Wall -fopenmp -c mst.c -o mst.o <br />

<br />

Then create the library files for the sort and data_structure files:

ar rcs libsort.a sort.o <br />
ar rcs libdata_structures.a data_structures.o <br />

<br />

Then, once we have compiled all the files, create the executable file:

mpicc -std=gnu99 -g -Wall -fopenmp -o mst mst.o data_structures.o sort.o -L<path_libraries> -lsort -ldata_structures


## Commands to run (inside an hpc cluster that works with jobs):

qsub run.sh

To run directly (not recommended to do this directly on the login node of the hpc cluster, but rather write this inside the .sh file)

mpirun.actual -np 4 <complete_path/mst> <complete_path/datasets/input_file.txt> 4
