#ifndef DATA_STRUCTURES
#define DATA_STRUCTURES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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
void newSet(Set* set, const int elements);

//initialize and allocate memory for the graph
void newGraph(Graph* graph, const int vertices, const int edges);

//return the parent of the input vertex using path compression
int findSet(Set* set, const int vertex);


//creates the Union of two sets with union by rank
void unionSet(Set* set, const int x, const int y);
#endif