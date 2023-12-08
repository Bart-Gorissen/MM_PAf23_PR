#ifndef PAGERANK_PARALLEL_H
#define PAGERANK_PARALLEL_H

#include "utility.h"
#include "CRSgraph.h"
#include "utility_parallel.h"

/**
 * Initializes random seed and parallel info
 * Return: PARInfo on BSP setup and data distribution
*/
PARInfo pagerank_par_init(long N, long seed);

/**
 * Parallel pagerank implementation
 * NOTE: naive: communicates full x to all processors during pGD^{-1}x
 * Return: local part of pagerank solution
*/
double* pagerank_par_naive(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V);

/**
 * Parallel pagerank implementation
 * NOTE: gget: communicates only elements needed to processors during pGD^{-1}x, but possibly with multiplicity
 * Return: local part of pagerank solution
*/
double* pagerank_par_gget(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V);

#endif /* PAGERANK_PARALLEL_H */