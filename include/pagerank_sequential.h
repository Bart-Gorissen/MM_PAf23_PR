#ifndef PAGERANK_SEQUENTIAL_H
#define PAGERANK_SEQUENTIAL_H

#include <stdio.h>
#include <math.h>

#include "utility.h"
#include "CRSgraph.h"

/**
 * Computes pagerank solution
 * pre_invert_D {0,1} : computes D^{-1} in memory
 * full_res {0,1} : ( 0 : computes r update from r ), ( 1 : computes r from u )
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq(CRSGraph graph, double p, double eps, int pre_invert_D, int full_res, int V);

#endif /* PAGERANK_SEQUENTIAL_H */