#ifndef CRSGRAPH_H
#define CRSGRAPH_H

#include <stdio.h>

// Characters used for printing adjacency matrix
static char MATRIX_EMPTY = '0';
static char MATRIX_FULL = '1';

// Graph datastructure (CRS representation)
typedef struct {
    long N; // # columns
    long M; // # rows (M <= N)
    long L; // max inlink #
    long* rowsize; // how many inlinks per node
    long* colindex; // inlink index (partitioned in rowsize)
    long nr_entries; // total number of inlinks added
} CRSGraph;

/**
 * Generates random graph on N nodes in CRS format
 * First determining how many entries per row: in range [1,L]
 * Then determining the location of each entry
 * 
 * Returns: fully set CRSGraph struct
*/
CRSGraph CRSGraph_generate(long N, long M, long L);


/**
 * Prints matrix representation of G
 * Return: void
*/
void CRSGraph_print(CRSGraph graph);


/**
 * Sorts colindex per row
 * Return: void
*/
void CRSGraph_sort(CRSGraph graph);

/**
 * Computes the column sum of the local adjacency matrix
 * Return: column sum
*/
long* CRSGraph_colsum(CRSGraph graph);

/**
 * Computes the column sum of the local adjacency matrix
 * NOTE: assumes colsum contains N zeros
 * Return: void
*/
void CRSGraph_colsum_inplace(CRSGraph graph, long* colsum);

/**
 * Compute stochstic-making diagonal
 * Return: column sum vector with 0 entries changed to 1 entries
*/
long* CRSGraph_stochastic_diagonal(CRSGraph graph);

/**
 * Computes the start indices of each row in the colindex array
 * Return: index list
*/
long* CRSGraph_indexlist(CRSGraph graph);

#endif /* CRSGRAPH_H */