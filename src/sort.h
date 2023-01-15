#ifndef SORT
#define SORT

//C header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <time.h>

//mpi and omp header files
#include <mpi.h>
#include <omp.h>

#include "sort.h"
#include "data_structures.h"

//merge sorted lists, start and end are inclusive
void merge(int* edgeList, const int start, const int end, const int pivot);

//sort the edge list using merge sort algorithm, start and end are inclusive
void mergeSort(int* edgeList, const int start, const int end);

//send the edge list of the graph from one process to all other processes in a communicator
void scatterEdgeList(int* edgeList, int* edgeListPart, const int elements,int* elementsPart);

// sort the edges of the graph in parallel with mergesort in parallel
void sort(Graph* graph);

#endif