#ifndef TESTS_H
#define TESTS_H

#include <stdlib.h>
#include <time.h>

#include "utility.h"
#include "CRSgraph.h"
#include "utility_parallel.h"
#include "pagerank_sequential.h"
#include "pagerank_parallel.h"

/**
 * Checks whether the random function produces values from [0, N)
 * This is done by repeated sampling, until all numbers are seen
 * NOTE: numbers of values outside the range are only reported on termination
 * Return: iters < max_iters
*/
int random_test(long N, long max_iters, int V);

/**
 * Tests graph construction, sorting, column sum, and stochastic diagonal
 * Return: 1
*/
int graph_seq_test(long N, long L, double p, int V);

/**
 * Prints the result of sorting CRS representation
 * Return: 1
*/
int graph_sort_seq_test(long N, long L, int V);

/**
 * Tests pagerank_seq function and scaling of resulting vector.
*/
int pagerank_seq_test1(long N, long L, double p, double eps, int V);

/**
 * Compares running time and number of iterations between sorted and unsorted CRS
*/
int pagerank_seq_test2(long N, long L, double p, double eps, int V);

/**
 * Compares running time and number of iterations between standard and pre-inverted D
*/
int pagerank_seq_test3(long N, long L, double p, double eps, int V);

/**
 * Compares iterations between fast and full residual computation
*/
int pagerank_seq_test4(long N, long L, double p, double eps, int V);

/**
 * Tests pagerank_par_naive function.
*/
int pagerank_par_test1(long N, long L, double p, double eps, int V);

/**
 * Tests pagerank_par_gget function.
*/
int pagerank_par_test2(long N, long L, double p, double eps, int V);

/**
 * Tests different computations for D function.
*/
int pagerank_par_test3(long N, long L, double p, double eps, int V);

#endif /* TESTS_H */