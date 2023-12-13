#include <math.h>
#include <bsp.h>

#include "utility_parallel.h"
#include "CRSgraph.h"

/**
 * Sets up PARInfo object for current processor
*/
PARInfo setup_parinfo(long N) {

    // set up the processor information
    PARInfo PI;
    PI.s = bsp_pid();
    PI.P = bsp_nprocs();

    // compute block size ceil( N / P )
    PI.b = N / PI.P;
    if (N % PI.P > 0) {
        PI.b ++;
    }

    if (N == 0) {
        PI.p_last = 0;
        PI.b_last = 0;
        PI.b_local = 0;
    }
    else {
        PI.p_last = N / PI.b; // number of full processors
        PI.b_last = N % PI.b; // elements to be put on first processor after full
        if (PI.s < PI.p_last) {
            PI.b_local = PI.b;
        }
        else if (PI.s > PI.p_last) {
            PI.b_local = 0;
        }
        else {
            PI.b_local = PI.b_last;
        }
    }

    return PI;
}

/**
 * Parallel sum of local_value on view PI, using par_vec to communicate
 * Uses: p data words (of type long) sent (and received), 1 synchronization
 * Return: sum_t local_value[@t]
*/
long parallel_sum_long(long local_value, long* par_vec, PARInfo PI) {

    // move the local value to each processor
    for (bsp_pid_t t=0; t<PI.P; t++) {
        bsp_put(t, &local_value, par_vec, PI.s * sizeof(long), sizeof(long));
    }

    bsp_sync();

    // add all local values
    long sum = 0;

    for (bsp_pid_t t=0; t<PI.P; t++) {
        sum += par_vec[t];
    }

    return sum;
}

/**
 * Parallel sum of local_value on view PI, using par_vec to communicate
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: sum_t local_value[@t]
*/
double parallel_sum(double local_value, double* par_vec, PARInfo PI) {

    // move the local value to each processor
    for (bsp_pid_t t=0; t<PI.P; t++) {
        bsp_put(t, &local_value, par_vec, PI.s * sizeof(double), sizeof(double));
    }

    bsp_sync();

    // add all local values
    double sum = 0;

    for (bsp_pid_t t=0; t<PI.P; t++) {
        sum += par_vec[t];
    }

    return sum;
}

/**
 * Parallel sum of squares (of vector)
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: sum_t sum_i (local_value[@t][i]^2)
*/
double parallel_sq_sum(double* local_vec, double* par_vec, PARInfo PI) {

    double local_value = 0;

    // first compute local sum of squares
    for (long i=0; i<PI.b_local; i++) {
        local_value += local_vec[i] * local_vec[i];
    }

    return parallel_sum(local_value, par_vec, PI);
}

/**
 * Parallel Euclidean norm
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: sqrt( sum_t sum_i (local_value[@t][i]^2) )
*/
double parallel_norm(double* local_vec, double* par_vec, PARInfo PI) {
    return sqrt(parallel_sq_sum(local_vec, par_vec, PI));
}

/**
 * Scales vector to stochastic
 * Uses: p data words (of type double) sent (and received), 1 synchronization
 * Return: global 1-norm prior to scaling
*/
double parallel_scale1(double* local_vec, double* par_vec, PARInfo PI) {
    double local_value = 0;
    for (long i=0; i<PI.b_local; i++) {
        local_value += local_vec[i];
    }

    double global_norm1 = parallel_sum(local_value, par_vec, PI);

    for (long i=0; i<PI.b_local; i++) {
        local_vec[i] /= global_norm1;
    }

    return global_norm1;
}

/**
 * Generates (stochastic) e_i vector
 * Return: e_i
*/
double* parallel_generate_ei(long i, PARInfo PI) {
    double* result = calloc(PI.b_local, sizeof(double));
    if (PI.b_local == 0) {
        return result;
    }
    long i_local = i - PI.s * PI.b;
    if (0 <= i_local && i_local < PI.b_local) {
        result[i] = 1;
    }
    return result;
}

/**
 * Column sum of graph using method D_comm_choice
 * D_comm_choice {0,1,2}: (0 : full-broadcast), (1 : P-round broadcast), (2 : message queue)
 * Return: colsum of PI.s part of graph
*/
long* parallel_colsum(CRSGraph graph, PARInfo PI, long* vec_par, int D_comm_choice) {
    long* colsum;

    switch (D_comm_choice) {
        case 0:
            colsum = parallel_colsum_bigpull(graph, PI, vec_par);
            break;
        
        case 1:
            // TODO
            colsum = parallel_colsum_bigpull(graph, PI, vec_par);
            break;

        case 2:
            // TODO
            colsum = parallel_colsum_bigpull(graph, PI, vec_par);
            break;
        
        default:
            colsum = NULL;
            break;
    }

    return colsum;
}

/**
 * Big-pull version of column sum of graph
 * NOTE: assumes vec_par has contains N zeros
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_bigpull(CRSGraph graph, PARInfo PI, long* vec_par) {
    CRSGraph_colsum_inplace(graph, vec_par);
    long* global_fetch = (long*) malloc(PI.b_local * PI.P * sizeof(long));

    long start_target = PI.s * PI.b * sizeof(long);
    long start_local = 0;

    for (bsp_pid_t t=0; t<PI.p_last; t++) {
        bsp_get(t,
                vec_par,
                start_target,
                &global_fetch[start_local],
                PI.b_local * sizeof(long)
        );

        start_local += PI.b_local;
    }

    bsp_sync();

    long* colsum = (long*) calloc(PI.b_local, sizeof(long));
    start_local = 0;

    for (bsp_pid_t t=0; t<PI.p_last; t++) {
        for (long i=0; i<PI.b_local; i++) {
            colsum[i] += global_fetch[start_local + i];
        }

        start_local += PI.b_local;
    }

    free(global_fetch);

    return colsum;
}

/**
 * Computes y = pGx using method pGx_comm_choice
 * pGx_comm_choice {0,1,2}: (0 : full-broadcast of u), (1 : P-round broadcast), (2 : quick-get), (3 : mapped-get)
 * Return: void
*/
void parallel_pGx(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI, int pGx_comm_choice) {
    switch (pGx_comm_choice) {
        case 0:
            parallel_pGx_bigpull(graph, p, x_par, y_vec, PI);
            break;
        
        case 1:
            // TODO (placeholder)
            parallel_pGx_bigpull(graph, p, x_par, y_vec, PI);
            break;
        
        case 2:
            parallel_pGx_gget(graph, p, x_par, y_vec, PI);
            break;

        case 3:
            // TODO (placeholder)
            parallel_pGx_gget(graph, p, x_par, y_vec, PI);
            break;
        
        default:
            break;
    }
}

/**
 * Big-pull version of y = p G_s . x
 * Return: void
*/
void parallel_pGx_bigpull(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI) {
    double* global_fetch = (double*) malloc(graph.N * sizeof(double));
    long start_local = 0;

    for (bsp_pid_t t=0; t<PI.p_last; t++) {
        bsp_get(t,
                x_par,
                0,
                &global_fetch[start_local],
                (t == PI.p_last) ? PI.b_last * sizeof(long) : PI.b_local * sizeof(long)
        );

        start_local += PI.b_local;
    }

    bsp_sync();

    long start = 0;

    for (long i=0; i<PI.b_local; i++) {
        y_vec[i] = 0;
        for (long j=0; j<graph.rowsize[i]; j++) {
            y_vec[i] += global_fetch[graph.colindex[start + j]];
        }

        y_vec[i] *= p;
        start += graph.rowsize[i];
    }

    free(global_fetch);
}

/**
 * G_get version of y = p G_s . x
 * Return: void
*/
void parallel_pGx_gget(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI) {
    double* global_fetch = (double*) malloc(graph.nr_entries * sizeof(double));

    // x_i is on proc t with elts t * b <= i < t * b_t_local
    // t = i / b
    // local index: i - t*b

    long start = 0;

    for (long i=0; i<PI.b_local; i++) {
        for (long j=0; j<graph.rowsize[i]; j++) {
            long k = graph.colindex[start + j];
            long t = k / PI.b;
            long k_local = k % PI.b;

            bsp_get(t,
                    x_par,
                    k_local * sizeof(double),
                    &global_fetch[start + j],
                    sizeof(double)
            );
        }

        start += graph.rowsize[i];
    }

    bsp_sync();

    start = 0;

    for (long i=0; i<PI.b_local; i++) {
        y_vec[i] = 0;
        for (long j=0; j<graph.rowsize[i]; j++) {
            y_vec[i] += global_fetch[start + j];
        }

        y_vec[i] *= p;
        start += graph.rowsize[i];
    }

    free(global_fetch);
}

/**
 * For init: we can save 1 synchronization, by combining generation of u_s and computation of D_s by splitting before and after sync part (as they are independent)
*/