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

typedef struct Set {
	int elements;
	int* parent;
	int* rank;
} Set;

typedef struct Graph {
	int numEdges;
	int numVertices;
	int* edgeList;
} Graph;

//initialize and allocate memory for the set
void newSet(Set* set, const int elements) {
	set->elements = elements;
	set->parent = (int*) malloc(elements * sizeof(int));
	memset(set->parent, -1, elements * sizeof(int));
	set->rank = (int*) calloc(elements, sizeof(int));
}

//initialize and allocate memory for the graph
void newGraph(Graph* graph, const int vertices, const int edges) {
	graph->numEdges = edges;
	graph->numVertices = vertices;
	graph->edgeList = (int*) calloc(edges * 3, sizeof(int));
}

//read input file to generate the input graph
void readInputFile(Graph* graph, const char inputFileName[]) {

	//open the file
	FILE* inputFile;
	inputFile = fopen(inputFileName, "r");
	if (inputFile == NULL) {
		fprintf(stderr, "Error: cannot open input file!\n");
		exit(1);
	}

	int fscanfResult;

	//first line contains the number of vertices and edges
	int vertices = 0;
	int edges = 0;
	fscanfResult = fscanf(inputFile, "%d %d", &vertices, &edges);

	//generate new graph with given number of vertices and edges
	newGraph(graph, vertices, edges);

	int from;
	int to;
	int weight;

	//generate edges
	for (int i = 0; i < edges; i++) {
		fscanfResult = fscanf(inputFile, "%d %d %d", &from, &to, &weight);
		graph->edgeList[i * 3] = from;
		graph->edgeList[i * 3 + 1] = to;
		graph->edgeList[i * 3 + 2] = weight;
	}

	fclose(inputFile);
}

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



//return the parent of the input vertex using path compression
int findSet(Set* set, int vertex) {
	if (set->parent[vertex] == -1) {
		return vertex;
	} 
	else {
		set->parent[vertex] = findSet(set,set->parent[vertex]);
		return set->parent[vertex];
	}
}


//creates the Union of two sets with union by rank
void unionSet(Set* set, const int x, const int y) {
	int xroot = findSet(set, x);
	int yroot = findSet(set, y);

	//first check if the two roots are the same
	if (xroot == yroot) {
		return;
	}
	//connecting tree with lowest rank to the tree with highest rank 
	else if (set->rank[xroot] < set->rank[yroot]) {
		set->parent[xroot] = yroot;
	} 
	else if (set->rank[xroot] > set->rank[yroot]) {
		set->parent[yroot] = xroot;
	} 
	//if ranks are same, arbitrarily increase the rank of one node
	else {
		set->parent[xroot] = yroot;
		set->rank[yroot] = set->rank[xroot] + 1;
	}
}



//find a MST of the graph using Kruskal's algorithm
void mstKruskal(Graph* graph, Graph* mst) {

	//create needed data structure
	Set* set = &(Set ) { .elements = 0, .parent = NULL, .rank =NULL };
	newSet(set, graph->numVertices);

	//get processor rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// sorting step
	sort(graph);

	//only one processor can execute the following block of code.
	//This because now the algorithm is adding the edges to the mst graph
	//and it is known that it is not possible to parallelize this part of
	//the algorithm since edges must be added one per time and it is
	//necessary to pay attention whenever adding an edge creates a cycle
	//or if two vertices are in the same subtree or not
	if (rank == 0) {

		//add edges to the MST; graph edge is the index of the edge 
		//in the edge list that is currently added to the graph.
		int graphEdge = 0;

		//for loop cycle; add edges until either all vertices have been added or
		//all edges have been added to the mst
		for (int mstEdge = 0; mstEdge < graph->numVertices - 1 || graphEdge < graph->numEdges;) {

			//check for loops if edge would be inserted
			int parentFrom = findSet(set, graph->edgeList[graphEdge * 3]);
			int parentTo = findSet(set, graph->edgeList[graphEdge * 3 + 1]);
			if (parentFrom != parentTo) {
				//add edge to MST
				memcpy(&mst->edgeList[mstEdge * 3], &graph->edgeList[graphEdge * 3], 3 * sizeof(int));
				unionSet(set, parentFrom, parentTo);
				mstEdge++;
			}
			graphEdge++;
		}
	}

	//cleanup
	free(set->parent);
	free(set->rank);
}



// main program
int main(int argc, char* argv[]) {

	// MPI initialization
	MPI_Init(&argc, &argv);

    //Get the number of processes
    int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Get the rank of the process
    int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	//creating graph variables. It is necessary to define two graphs:
    //1) one for the input graph
    //2) one for the mst (which is a sub-graph of the input graph)
	Graph* graph = &(Graph) { .numEdges = 0, .numVertices = 0, .edgeList = NULL };
	Graph* mst = &(Graph) { .numEdges = 0, .numVertices = 0, .edgeList = NULL };

    //only one processor has to read the input graph
	if (rank == 0) {

		//read the file and generate input graph
		readInputFile(graph, argv[1]);

		//generate the mst graph
		newGraph(mst, graph->numVertices, graph->numVertices - 1);
	}

	//Create the start to compute how much time it takes to find the mst of the given graph
	//using Kruskal's algorithm  
	double start_time = MPI_Wtime();
	mstKruskal(graph, mst);
	double finish_time;

	//only one processor has to print the results
	if (rank == 0) {

		//compute the mst weight
		//here it is possible to use omp
		unsigned long mstWeight = 0;
		#pragma omp parallel for reduction(+: mstWeight)
		for (int i = 0; i < mst->numEdges; i++) {
			mstWeight += mst->edgeList[i * 3 + 2];
		}

		//define the interval time took to find the mst
		finish_time = MPI_Wtime();

		//print the edges of the MST, its weight and the computation time
		printf("Minimum Spanning Tree:\n");
		printf("------------------------------------------------\n");
		for (int i = 0; i < mst->numEdges; i++) {
			for (int j = 0; j < 3; j++) {
				printf("%d\t", mst->edgeList[i * 3 + j]);
			}
			printf("\n");
		}
		printf("------------------------------------------------\n");
		printf("MST weight: %lu\n", mstWeight);
		printf("Time (sec): %f s\n", finish_time - start_time);

		//cleanup the two graphs
		free(graph->edgeList);
		free(mst->edgeList);
	}

    // Finalize the MPI environment
	MPI_Finalize();

	return 0;
}