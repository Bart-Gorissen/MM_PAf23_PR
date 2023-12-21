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
 * Column sum of graph using method D_comm_choice
 * D_comm_choice {0,1,2,3}: (0 : full-broadcast), (1 : P-round broadcast), (2 : message queue), (3 : message queue - aggregated)
 * Return: colsum of PI.s part of graph
*/
long* parallel_colsum(CRSGraph graph, PARInfo PI, long* vec_par, int D_comm_choice);

/**
 * Big-pull version of column sum of graph
 * Return: colsum for PI.s part of graph
*/
long* parallel_colsum_bigpull(CRSGraph graph, PARInfo PI, long* vec_par);

/**
 * P-Round version of column sum of graph
 * NOTE: Assumes vec_par has PI.b entries
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_prounds(CRSGraph graph, PARInfo PI, long* vec_par);

/**
 * Blind-send version of column sum of graph
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_blindsend(CRSGraph graph, PARInfo PI);

/**
 * Aggregated-send version of column sum of graph
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_aggrsend(CRSGraph graph, PARInfo PI);

/**
 * Computes y = pGx using method pGx_comm_choice
 * pGx_comm_choice {0,1,2,3}: (0 : full-broadcast of u), (1 : P-round broadcast), (2 : quick-get), (3 : mapped-get)
 * Return: void
*/
void parallel_pGx(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI, int pGx_comm_choice, IndexMap IM);

/**
 * Big-pull version of y = p G_s . x
 * Return: void
*/
void parallel_pGx_bigpull(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI);

/**
 * P-Round version of y = p G_s . x
 * NOTE: assumes graph is already sorted
 * Return: void
*/
long* parallel_pGx_pround(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI);

/**
 * G_get version of y = p G_s . x
 * Return: void
*/
void parallel_pGx_gget(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI);

/**
 * Mapped Get version of y = p G_s . x
 * Return: void
*/
void parallel_pGx_mapget(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI, IndexMap IM);

#endif /* UTILITY_PARALLEL_H */