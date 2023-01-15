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

//include the other files
#include "sort.h"
#include "data_structures.h"

//merge sorted lists, start and end are inclusive
void merge(int* edgeList, const int start, const int end, const int pivot) {
	int length = end - start + 1;
	int* working = (int*) malloc(length * 3 * sizeof(int));

	//copy first part
	memcpy(working, &edgeList[start * 3],(pivot - start + 1) * 3 * sizeof(int));

	//copy second part
	int workingEnd = end + pivot - start + 1;
	for (int i = pivot + 1; i <= end; i++) {
		memcpy(&working[(workingEnd - i) * 3], &edgeList[i * 3], 3 * sizeof(int));
	}

	int left = 0;
	int right = end - start;
	for (int k = start; k <= end; k++) {
		if (working[right * 3 + 2]< working[left * 3 + 2]) {
			memcpy(&edgeList[k * 3], &working[right * 3], 3 * sizeof(int));
			right--;
		} else {
			memcpy(&edgeList[k * 3], &working[left * 3], 3 * sizeof(int));
			left++;
		}
	}

	//clean up
	free(working);
}

//sort the edge list using merge sort algorithm, start and end are inclusive
void mergeSort(int* edgeList, const int start, const int end) {
	if (start != end) {
		//recursively divide the list in two parts and sort them
		int pivot = (start + end) / 2;
		mergeSort(edgeList, start, pivot);
		mergeSort(edgeList, pivot + 1, end);

		merge(edgeList, start, end, pivot);
	}
}

//send the edge list of the graph from one process to all other processes in a communicator
void scatterEdgeList(int* edgeList, int* edgeListPart, const int elements,int* elementsPart) {

    //Get the number of processes
    int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Get the rank of the process
    int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Scatter(edgeList, *elementsPart * 3, MPI_INT, edgeListPart,	*elementsPart * 3, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == size - 1 && elements % *elementsPart != 0) {
		//number of elements and processes isn't divisible without remainder
		//it is necessary to do this only for one process
		*elementsPart = elements % *elementsPart;
	}

	if (elements / 2 + 1 < size && elements != size) {
		if (rank == 0) {
			fprintf(stderr, "Unsupported size/process combination, exiting!\n");
		}
		MPI_Finalize();
		exit(1);
	}
}

// sort the edges of the graph in parallel with mergesort in parallel
void sort(Graph* graph) {

    //Get the number of processes
    int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Get the rank of the process
    int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	//send number of elements to the other processes 
	int elements;
	if (rank == 0) {
		elements = graph->numEdges;
		MPI_Bcast(&elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} 
	else {
		MPI_Bcast(&elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	//Sends the edges to sort from one process to all other processes in a communicator
	int elementsPart = (elements + size - 1) / size;
	int* edgeListPart = (int*) malloc(elementsPart * 3 * sizeof(int));
	scatterEdgeList(graph->edgeList, edgeListPart, elements, &elementsPart);

	//sorting step 
	mergeSort(edgeListPart, 0, elementsPart - 1);

	//merge all parts
	int from;
	int to;
	int elementsRecieved;
	for (int step = 1; step < size; step *= 2) {
		if (rank % (2 * step) == 0) {
			from = rank + step;
			if (from < size) {
				MPI_Recv(&elementsRecieved, 1, MPI_INT, from, 0,MPI_COMM_WORLD, &status);
				edgeListPart = realloc(edgeListPart,(elementsPart + elementsRecieved) * 3* sizeof(int));
				MPI_Recv(&edgeListPart[elementsPart * 3],elementsRecieved * 3,MPI_INT, from, 0, MPI_COMM_WORLD, &status);
				merge(edgeListPart, 0, elementsPart + elementsRecieved - 1,	elementsPart - 1);
				elementsPart += elementsRecieved;
			}
		} 
		else if (rank % step == 0) {
			to = rank - step;
			MPI_Send(&elementsPart, 1, MPI_INT, to, 0, MPI_COMM_WORLD);
			MPI_Send(edgeListPart, elementsPart * 3, MPI_INT, to,0,	MPI_COMM_WORLD);
		}
	}

	//edgeListPart is the new edgeList of the graph 
	//cleanup other memory
	if (rank == 0) {
		free(graph->edgeList);
		graph->edgeList = edgeListPart;
	} else {
		free(edgeListPart);
	}
}