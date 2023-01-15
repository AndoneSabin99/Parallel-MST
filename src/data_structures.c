#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "data_structures.h"

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

//return the parent of the input vertex using path compression
int findSet(Set* set, const int vertex) {
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