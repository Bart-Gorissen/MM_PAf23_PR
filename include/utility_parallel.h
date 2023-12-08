#ifndef UTILITY_PARALLEL_H
#define UTILITY_PARALLEL_H

#include <stdlib.h>
#include <math.h>
#include <bsp.h>

#include "utility_parallel.h"
#include "CRSgraph.h"

// Parallel execution information
typedef struct {
    bsp_pid_t s;        // processor number
    long b_local;       // processor s block size
    bsp_pid_t P;        // total number of processors
    long b;             // block size (general)
    bsp_pid_t p_last;   // last processor with non-zero block
    long b_last;        // block size (last processor with non-zero)
} PARInfo;

/**
 * Sets up PARInfo object for current processor
*/
PARInfo setup_parinfo(long N);

/**
 * Parallel sum of local_value on view PI, using par_vec to communicate
 * Uses: p data words (of type long) sent (and received), 1 synchronization
 * Return: sum_t local_value[@t]
*/
long parallel_sum_long(long local_value, long* par_vec, PARInfo PI);

/**
 * Parallel sum of local_value on view PI, using par_vec to communicate
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: sum_t local_value[@t]
*/
double parallel_sum(double local_value, double* par_vec, PARInfo PI);

/**
 * Parallel sum of squares (of vector)
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: sum_t sum_i (local_value[@t][i]^2)
*/
double parallel_sq_sum(double* local_vec, double* par_vec, PARInfo PI);

/**
 * Parallel Euclidean norm
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: sqrt( sum_t sum_i (local_value[@t][i]^2) )
*/
double parallel_norm(double* local_vec, double* par_vec, PARInfo PI);

/**
 * Scales vector to stochastic
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: global 1-norm prior to scaling
*/
double parallel_scale1(double* local_vec, double* par_vec, PARInfo PI);

/**
 * Generates (stochastic) e_i vector
 * Return: e_i
*/
double* parallel_generate_ei(long i, PARInfo PI);

/**
 * Big-pull version of column sum of graph
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_bigpull(CRSGraph graph, PARInfo PI);

/**
 * Big-pull version of y = p G_s . x
 * Return: void
*/
void parallel_pGy_bigpull(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI);

/**
 * G_get version of y = p G_s . x
 * Return: void
*/
void parallel_pGy_gget(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI);

#endif /* UTILITY_PARALLEL_H */