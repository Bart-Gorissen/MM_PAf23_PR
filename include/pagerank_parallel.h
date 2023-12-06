#ifndef PAGERANK_PARALLEL_H
#define PAGERANK_PARALLEL_H

#include "utility.h"
#include "CRSgraph.h"
#include "utility_parallel.h"

PARInfo pagerank_par_init(long N, long seed);

double* pagerank_par_naive(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V);

#endif /* PAGERANK_PARALLEL_H */