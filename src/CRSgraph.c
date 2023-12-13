#include <stdio.h>

#include "utility.h"
#include "CRSgraph.h"

/**
 * Generates random graph on N nodes in CRS format
 * First determining how many entries per row: in range [1,L]
 * Then determining the location of each entry
 * 
 * Returns: fully set CRSGraph struct
*/
CRSGraph CRSGraph_generate(long N, long M, long L) {
    CRSGraph graph;
    graph.N = N;
    graph.M = M;
    graph.L = L;

    // allocate the row size array
    graph.rowsize = (long*) malloc(N * sizeof(long));
    graph.nr_entries = 0;

    // determine row size
    for (long i=0; i<M; i++) {
        graph.rowsize[i] = 1 + urandom_rrange(L);
        graph.nr_entries += graph.rowsize[i];
    }

    // allocate space for the required number of entries
    graph.colindex = (long*) malloc(graph.nr_entries * sizeof(long));

    // for each row determine columns present
    for (long i=0; i<graph.nr_entries; i++) {
        graph.colindex[i] = urandom_rrange(N);
    }

    return graph;
}

/**
 * Prints matrix representation of G
 * Return: void
*/
void CRSGraph_print(CRSGraph graph) {
    // every row has N+1 characters, the last being "\n"
    char* result = (char*) malloc((graph.M*(graph.N+1) + 1) * sizeof(char));
    result[graph.M*(graph.N+1)] = '\0';

    long start_result = 0;
    long start_colindex = 0;

    // set all entries to [not-set]
    for(long i=0; i<graph.M*(graph.N+1); i++) {
        result[i] = MATRIX_EMPTY;
    }

    // fill-in matrix entries
    for (long i=0; i<graph.M; i++) {
        for (long j=0; j<graph.rowsize[i]; j++) {
            result[start_result + graph.colindex[start_colindex + j]] = MATRIX_FULL;            
        }
        result[start_result+graph.N] = '\n';

        start_result += graph.N+1;
        start_colindex += graph.rowsize[i];
    }

    // print graph information
    printf("Graph [%ld nodes, %ld rows, %ld max inlink, %ld entries]:\n%s", graph.N, graph.M, graph.L, graph.nr_entries, result);
    fflush(stdout);

    free(result);
}

/**
 * Sorts colindex per row
 * Return: void
*/
void CRSGraph_sort(CRSGraph graph) {
    if (graph.N < 1) return;

    long start = 0;

    for (long i=0; i<graph.M; i++) {
        quicksort(graph.colindex, start, start+graph.rowsize[i]-1);
        start += graph.rowsize[i];
    }
}

/**
 * Computes the column sum of the local adjacency matrix
 * Return: column sum
*/
long* CRSGraph_colsum(CRSGraph graph) {
    long* colsum = (long*) calloc(graph.N, sizeof(long));
    CRSGraph_colsum_inplace(graph, colsum);

    return colsum;
}

/**
 * Computes the column sum of the local adjacency matrix
 * NOTE: assumes colsum contains N zeros
 * Return: void
*/
void CRSGraph_colsum_inplace(CRSGraph graph, long* colsum) {
    long start = 0;

    for (long i=0; i<graph.N; i++) {
        colsum[i] = 0;
    }

    // add number of references to each column
    for (long i=0; i<graph.M; i++) {
        for (long j=0; j<graph.rowsize[i]; j++) {
            colsum[graph.colindex[start + j]] ++;
        }
        start += graph.rowsize[i];
    }
}

/**
 * Compute stochstic-making diagonal
 * Return: column sum vector with 0 entries changed to 1 entries
*/
long* CRSGraph_stochastic_diagonal(CRSGraph graph) {
    long* diagonal = CRSGraph_colsum(graph);
    change_01(diagonal, graph.N);

    return diagonal;
}

/**
 * Computes the start indices of each row in the colindex array (including one entry for the number of entries)
 * Return: index list
*/
long* CRSGraph_indexlist(CRSGraph graph) {
    long* indices = (long*) calloc(graph.M+1, sizeof(long));

    // set indices to first entry of each row
    for (long i=0; i<graph.M; i++) {
        for (long j=i+1; j<=graph.M; j++) {
            indices[j] += graph.rowsize[i];
        }
    }

    return indices;
}