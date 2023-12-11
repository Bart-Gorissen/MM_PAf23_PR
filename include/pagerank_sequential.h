#ifndef PAGERANK_SEQUENTIAL_H
#define PAGERANK_SEQUENTIAL_H

#include <stdio.h>
#include <math.h>

#include "utility.h"
#include "CRSgraph.h"

/**
 * Computes pagerank solution
 * NOTE: pre-inverts matrix D if pre_invert_D != 0
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq(CRSGraph graph, double p, double eps, int pre_invert_D, int V);

/**
 * Computes pagerank solution
 * NOTE: does not pre-invert matrix D
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq_initialD(CRSGraph graph, double p, double eps, int V);

/**
 * Computes pagerank solution
 * NOTE: pre-inverts matrix D
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq_invertedD(CRSGraph graph, double p, double eps, int V);

/**
 * Computes pagerank solution (using full computation for residual)
 * NOTE: pre-inverts matrix D if pre_invert_D != 0
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq_fullres(CRSGraph graph, double p, double eps, int pre_invert_D, int V);

/**
 * Computes pagerank solution (using full computation for residual)
 * NOTE: does not pre-invert matrix D
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq_initialD_fullres(CRSGraph graph, double p, double eps, int V);

/**
 * Computes pagerank solution (using full computation for residual)
 * NOTE: pre-inverts matrix D
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq_invertedD_fullres(CRSGraph graph, double p, double eps, int V);

#endif /* PAGERANK_SEQUENTIAL_H */