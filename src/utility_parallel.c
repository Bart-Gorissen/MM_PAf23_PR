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

    

    // PI.b_local = PI.b;
    // PI.b_last = N - (PI.b * (PI.P - 1));
    // if (PI.s == PI.P - 1) {
    //     PI.b_local = PI.b_last;
    // }

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
 * Big-pull version of column sum of graph
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_bigpull(CRSGraph graph, PARInfo PI) {
    long* local_vec = CRSGraph_colsum(graph);
    long* global_fetch = (long*) malloc(PI.b_local * PI.P * sizeof(long));

    bsp_push_reg(local_vec, graph.N * sizeof(long));
    bsp_sync();

    long start_target = PI.s * PI.b * sizeof(long);
    long start_local = 0;

    //printf("Proc %d, %ld, %ld, %ld\n", PI.s, PI.b_local, PI.b, PI.P); fflush(stdout);

    for (bsp_pid_t t=0; t<PI.p_last; t++) {
        bsp_get(t,
                local_vec,
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

    bsp_pop_reg(local_vec);
    free(local_vec);
    free(global_fetch);

    return colsum;
}

/**
 * Big-pull version of y = p G_s . x
 * Return: void
*/
void parallel_pGy_bigpull(CRSGraph graph, double p, double* x_par, double* y_vec, PARInfo PI) {
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
 * For init: we can save 1 synchronization, by combining generation of u_s and computation of D_s by splitting before and after sync part (as they are independent)
*/