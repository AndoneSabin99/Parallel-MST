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

	//free the set data structure
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

		//define the interval time took to find the mst and to compute its cost in terms of weight
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

		//free the two graphs
		free(graph->edgeList);
		free(mst->edgeList);
	}

    // Finalize the MPI environment
	MPI_Finalize();

	return 0;
}