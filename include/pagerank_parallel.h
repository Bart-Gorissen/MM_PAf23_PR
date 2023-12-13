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
 * pre_invert_D {0,1}: pre-computing the inverse of D
 * u_choice {0,1,...,N}: (0 : random), (i : e_i)
 * D_comm_choice {0,1,2}: (0 : full-broadcast), (1 : P-round broadcast), (2 : message queue)
 * pGy_comm_choice {0,1,2}: (0 : full-broadcast of u), (1 : P-round broadcast), (2 : quick-get), (3 : mapped-get)
 * Return: local part of pagerank solution
*/
double* pagerank_par(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int u_choice, int D_comm_choice, int pGy_comm_choice, int V);

// /**
//  * Parallel pagerank implementation
//  * NOTE: naive: communicates full x to all processors during pGD^{-1}x
//  * Return: local part of pagerank solution
// */
// double* pagerank_par_naive(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V);

// /**
//  * Parallel pagerank implementation
//  * NOTE: gget: communicates only elements needed to processors during pGD^{-1}x, but possibly with multiplicity
//  * Return: local part of pagerank solution
// */
// double* pagerank_par_gget(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V);

#endif /* PAGERANK_PARALLEL_H */