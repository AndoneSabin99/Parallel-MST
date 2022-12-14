#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>

typedef struct edge {
    int u, v, w;
} edge;

int n;
edge* edges, *chosen_edges;
int* par;
int num_edge;

int find_set(int u) {
    if (par[u] == u) return u;
    return par[u] = find_set(par[u]);
}

int merge_set(int u, int v) {
    int pu = find_set(u), pv = find_set(v);
    if (pu == pv) return 0;
    par[pv] = pu;
    return 1;
}

int comparison_weight(edge* x, edge* y) {
    if (x->w == y->w) {
        if (x->u == y->u)
            return x->v < y->v;
        return x->u < y->u;
    }
    return x->w < y->w;
}

int comparison_node(edge* x, edge* y) {
    if (x->u == y->u)
        return x->v < y->v;
    return x->u < y->u;
}

void merge(edge edges[], edge larr[], int nl, edge rarr[], int nr, int (*comparison)(edge*, edge*), int offset) {
    int il = 0, ir = 0, j = offset;
    while (il < nl && ir < nr) {
        if ((*comparison)(&larr[il], &rarr[ir])) {
            edges[j] = larr[il];
            il++;
        } else {
            edges[j] = rarr[ir];
            ir++;
        }
        j++;
    }

    while (il < nl) {
        edges[j] = larr[il];
        il++; j++;
    }

    while (ir < nr) {
        edges[j] = rarr[ir];
        ir++; j++;
    }
}

void merge_sort(edge edges[], int n, int (*comparison)(edge*, edge*)) {
    if (n > 1) {
        int m = n / 2;
        edge* larr = malloc(m * sizeof(edge));
        edge* rarr = malloc((n - m) * sizeof(edge));
        memcpy(larr, edges, m * sizeof(edge));
        memcpy(rarr, edges + m, (n - m) * sizeof(edge));

        #pragma omp parallel
        {
            #pragma omp single
            {
                #pragma omp task shared(larr, m, comparison)
                merge_sort(larr, m, comparison);

                #pragma omp task shared(rarr, n, m, comparison)
                merge_sort(rarr, n - m, comparison);
            }
        }

        merge(edges, larr, m, rarr, n - m, comparison, 0);
        free(larr); free(rarr);
    }
}

void merge_gather(edge gather[], int sendcounts[], int displ[], int (*comparison)(edge*, edge*), int l, int r) {
    if (l >= r) return;
    int mid = (l + r) >> 1;
    //TO DO: need to check if we can apply openmp also here
    merge_gather(gather, sendcounts, displ, comparison, l, mid);
    merge_gather(gather, sendcounts, displ, comparison, mid + 1, r);
    int nl = 0, nr = 0;
    for (int i = l; i <= r; i++) {
        if (i <= mid)
            nl += sendcounts[i];
        else
            nr += sendcounts[i];
    }
    edge* larr = malloc(nl * sizeof(edge));
    edge* rarr = malloc(nr * sizeof(edge));
    memcpy(larr, gather + displ[l], nl * sizeof(edge));
    memcpy(rarr, gather + displ[l] + nl, nr * sizeof(edge));
    merge(gather, larr, nl, rarr, nr, comparison, displ[l]);
    free(larr); free(rarr);
}

int main(int argc, char** argv) {
    int world_size, world_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Define mpi type
    MPI_Datatype MPI_EDGE;
    // Type create struct
    int block_count = 3;
    const int block_lengths[3] = { 1, 1, 1 };
    const MPI_Aint displacements[3] = { 0, sizeof(int), sizeof(int) * 2 };
    MPI_Datatype block_types[3] = { MPI_INT, MPI_INT, MPI_INT };
    MPI_Type_create_struct(block_count, block_lengths, displacements, block_types, &MPI_EDGE);
    // Commit type
    MPI_Type_commit(&MPI_EDGE);


    //getting number of threads for OpenMP
    int num_thread = omp_get_max_threads();
    // omp_set_max_active_levels(omp_get_max_threads());

    // Input
    //TO DO : need to change this according to our project, either with .csv or .txt file
    //in this case just read from a file the number of nodes and the values of weights

    FILE *graphFile;
    graphFile = fopen(argv[1], "r");

    if (graphFile == NULL){
        printf("Error Reading File\n");
        exit (0);
    }

    //t variable is used in order to track the time from where the program starts its execution
    clock_t t = clock();
    if (world_rank == 0) {
        int n, x;
        fscanf(graphFile, "%d", &n);
        edges = (edge*) malloc(n * (n + 1) / 2 * sizeof(edge));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int x;
                fscanf(graphFile, "%d", &x);
                if (x == -1) continue;
                if (i >= j) continue;
                edges[num_edge].u = i;
                edges[num_edge].v = j;
                edges[num_edge].w = x;
                num_edge++;
            }
        }
        assert(num_edge >= n - 1);
    }
    fclose(graphFile);
    MPI_Barrier(MPI_COMM_WORLD);

    //send the input to the other nodes
    if (world_rank == 0) {
        for (int i = 1; i < world_size; i++)
            MPI_Send(&num_edge, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&num_edge, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    

    // Process
    {
        // Define array for gather and scatter
        edge* gathered = malloc(num_edge * sizeof(edge));
        int scattered_size = num_edge / world_size;
        edge* scattered = malloc((scattered_size + 1) * sizeof(edge));

        int* sendcounts = malloc(world_size * sizeof(int));
        int* displs = malloc(world_size * sizeof(int));
        int rem = num_edge % world_size;
        int sum = 0;
        for (int i = 0; i < world_size; i++) {
            sendcounts[i] = num_edge / world_size + (i < rem);
            displs[i] = sum;
            sum += sendcounts[i];
        }

        // Scatter and gather
        MPI_Scatterv(edges, sendcounts, displs, MPI_EDGE, scattered, sendcounts[world_rank], MPI_EDGE, 0, MPI_COMM_WORLD);
        merge_sort(scattered, sendcounts[world_rank], comparison_weight);
        MPI_Gatherv(scattered, sendcounts[world_rank], MPI_EDGE, gathered, sendcounts, displs, MPI_EDGE, 0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            for (int i = 0; i < num_edge; i++) {
                edges[i] = gathered[i];
            }
            merge_gather(gathered, sendcounts, displs, comparison_weight, 0, world_size - 1);
            for (int i = 0; i < num_edge; i++) {
                edges[i] = gathered[i];
            }
        }
        free(gathered); free(scattered);
        free(sendcounts); free(displs);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int num_chosen = 0;
    long long total_cost = 0;
    if (world_rank == 0) {
        par = (int*) malloc(n * sizeof(int));
        for (int i = 0; i < n; i++) {
            par[i] = i;
        }
        chosen_edges = (edge*) malloc(num_edge * sizeof(edge));
        for (int i = 0; i < num_edge; i++) {
            int u = edges[i].u;
            int v = edges[i].v;
            int w = edges[i].w;
            if (merge_set(u, v)) {
                total_cost += w;
                chosen_edges[num_chosen++] = edges[i];
                if (num_chosen == n - 1) break;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        for (int i = 1; i < world_size; i++)
            MPI_Send(&num_chosen, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&num_chosen, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    {
        // Define array for gather and scatter
        edge* gathered = malloc(num_chosen * sizeof(edge));
        int scattered_size = num_chosen / world_size;
        edge* scattered = malloc((scattered_size + 1) * sizeof(edge));

        int* sendcounts = malloc(world_size * sizeof(int));
        int* displs = malloc(world_size * sizeof(int));
        int rem = num_chosen % world_size;
        int sum = 0;
        for (int i = 0; i < world_size; i++) {
            sendcounts[i] = num_chosen / world_size + (i < rem);
            displs[i] = sum;
            sum += sendcounts[i];
        }

        // Scatter and gather
        MPI_Scatterv(chosen_edges, sendcounts, displs, MPI_EDGE, scattered, sendcounts[world_rank], MPI_EDGE, 0, MPI_COMM_WORLD);
        merge_sort(scattered, sendcounts[world_rank], comparison_node);
        MPI_Gatherv(scattered, sendcounts[world_rank], MPI_EDGE, gathered, sendcounts, displs, MPI_EDGE, 0, MPI_COMM_WORLD);

         if (world_rank == 0) {
            merge_gather(gathered, sendcounts, displs, comparison_node, 0, world_size - 1);
            for (int i = 0; i < num_chosen; i++) {
                chosen_edges[i] = gathered[i];
            }
        }

        free(gathered); free(scattered);
        free(sendcounts); free(displs);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Output
    if (world_rank == 0) {
        printf("%lld\n", total_cost);

        for (int i = 0; i < num_chosen; i++) {
            printf("%d-%d\n", chosen_edges[i].u, chosen_edges[i].v);
        }
        double time_taken = ((double) (clock() - t)) / CLOCKS_PER_SEC;
        printf("Time taken: %f ms\n", time_taken);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // End
    MPI_Type_free(&MPI_EDGE);
    MPI_Finalize();
    return 0;
}