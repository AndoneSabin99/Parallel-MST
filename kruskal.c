#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// main program
int main(int argc, char* argv[]) {
	// MPI variables and initialization
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// initialize the graph

	if (rank == 0) {
    
    //finish the initialization of the graph here

		// (eventually) print the edges of the graph to be sure it is ok

		
	}

  
	double start = MPI_Wtime();
	// use Kruskal's algorithm
	//mstKruskal(graph, mst);
	

	if (rank == 0) {

		// print the edges of the MST
		printf("Minimum Spanning Tree (Kruskal):\n");
		//printWeightedGraph(mst);

    /*
    //compute weightMST
		unsigned long weightMST = 0;
		for (int i = 0; i < mst->edges; i++) {
			weightMST += mst->edgeList[i * 3 + 2];
		}
		printf("MST weight: %lu\n", weightMST);
    */

		// cleanup the graph

		printf("Time elapsed: %f s\n", MPI_Wtime() - start);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
