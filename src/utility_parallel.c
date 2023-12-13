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

    if (N % PI.P == 0) {
        PI.b = N / PI.P;
        PI.p_last = PI.P-1;
        PI.b_local = PI.b;
        PI.b_last = PI.b;
    }
    else {
        PI.b = (N / PI.P) + 1;
        PI.p_last = N / PI.b;
        PI.b_last = N % PI.b;
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
            colsum = parallel_colsum_prounds(graph, PI, vec_par);
            break;

        case 2:
            colsum = parallel_colsum_blindsend(graph, PI);
            break;

        case 3:
            colsum = parallel_colsum_aggrsend(graph, PI);
            break;

        default:
            colsum = NULL;
            break;
    }

    return colsum;
}

/**
 * Big-pull version of column sum of graph
 * NOTE: assumes vec_par contains N zeros
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_bigpull(CRSGraph graph, PARInfo PI, long* vec_par) {
    CRSGraph_colsum_inplace(graph, vec_par);
    long* global_fetch = (long*) malloc(PI.b_local * (PI.p_last+1) * sizeof(long));

    long start_target = PI.s * PI.b * sizeof(long);
    long start_local = 0;

    for (bsp_pid_t t=0; t<=PI.p_last; t++) {
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

    for (bsp_pid_t t=0; t<=PI.p_last; t++) {
        for (long i=0; i<PI.b_local; i++) {
            colsum[i] += global_fetch[start_local + i];
        }

        start_local += PI.b_local;
    }

    free(global_fetch);

    return colsum;
}

/**
 * P-Round versoin of column sum of graph
 * NOTE: Assumes vec_par conts PI.b zeros
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_prounds(CRSGraph graph, PARInfo PI, long* vec_par) {
    long* colsum = (long*) calloc(PI.b_local, sizeof(long));
    long* global_fetch = (long*) malloc(PI.b_local * sizeof(long));

    // method assumes sorted graph
    CRSGraph_sort(graph);

    long* start_indices = CRSGraph_indexlist(graph);
    long* curr_indices = (long*) malloc(PI.b_local * sizeof(long));

    // copy starting indices and forward to start current processor
    for (long i=0; i<PI.b_local; i++) {
        curr_indices[i] = start_indices[i];

        // forward
        while(graph.colindex[curr_indices[i]] / PI.b < PI.s && curr_indices[i] < start_indices[i+1]) {
            curr_indices[i] ++;
        }
        vec_par[i] = 0;
    }
    for (long i=PI.b_local; i<PI.b; i++) {
        vec_par[i] = 0;
    }

    for (bsp_pid_t t=0; t<=PI.p_last; t++) {
        
        bsp_pid_t target = (PI.s + t) % (PI.p_last + 1);
        bsp_pid_t source = (PI.s + (PI.p_last + 1) - t) % (PI.p_last + 1);

        // wrap around when applicable (not needed for processor 0)
        if (target == 0 && PI.s != 0) {
            for (long i=0; i<PI.b_local; i++) {
                curr_indices[i] = start_indices[i];
            }
        }

        for (long i=0; i<PI.b_local; i++) {
            // while the next column is on the target processor, increment the appropriate entry
            while (graph.colindex[curr_indices[i]] / PI.b <= target && curr_indices[i] < start_indices[i+1]) {
                vec_par[graph.colindex[curr_indices[i]] % PI.b] ++;
                curr_indices[i] ++;
            }
        }

        // get the sum ()
        bsp_get(source,
                vec_par,
                0,
                global_fetch,
                PI.b_local * sizeof(long)
        );

        bsp_sync();

        // increment the sum and set the (running) columns sum to 0
        for (long i=0; i<PI.b_local; i++) {
            colsum[i] += global_fetch[i];
            vec_par[i] = 0;
        }
        for (long i=PI.b_local; i<PI.b; i++) {
            vec_par[i] = 0;
        }
    }

    free(global_fetch);
    free(start_indices);
    free(curr_indices);

    return colsum;
}

/**
 * Blind-send version of column sum of graph
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_blindsend(CRSGraph graph, PARInfo PI) {
    long* colsum = (long*) calloc(PI.b_local, sizeof(long));

    long start = 0;

    // send all columns indices to their responsible processors
    for (long i=0; i<PI.b_local; i++) {
        for (long j=0; j<graph.rowsize[i]; j++) {
            bsp_send(graph.colindex[start + j] / PI.b,
                     NULL,
                     &graph.colindex[start + j],
                     sizeof(long)
            );
        }

        start += graph.rowsize[i];
    }

    bsp_sync();

    // determine number of messages received
    int nr_msgs;
    bsp_qsize(&nr_msgs, NULL);

    // for each message, increment the column
    for (long m=0; m<nr_msgs; m++) {
        long col_global;
        bsp_move(&col_global, sizeof(long));
        colsum[col_global % PI.b] ++;
    }

    return colsum;
}

/**
 * Aggregated-send version of column sum of graph
 * Return: colsum for PI.s of graph
*/
long* parallel_colsum_aggrsend(CRSGraph graph, PARInfo PI) {
    long* colsum = (long*) calloc(PI.b_local, sizeof(long));

    // method assumes sorted graph
    CRSGraph_sort(graph);

    long* start_indices = CRSGraph_indexlist(graph);
    long* curr_indices = (long*) calloc(PI.b, sizeof(long));
    long* curr_colsum = (long*) calloc(PI.b, sizeof(long));

    // copy starting indices
    for (long i=0; i<PI.b_local; i++) {
        curr_indices[i] = start_indices[i];
    }

    for (bsp_pid_t t=0; t<=PI.p_last; t++) {
        for (long i=0; i<PI.b_local; i++) {
            while (graph.colindex[curr_indices[i]] / PI.b <= t && curr_indices[i] < start_indices[i+1]) {
                curr_colsum[graph.colindex[curr_indices[i]] % PI.b] ++;
                curr_indices[i] ++;
            }
        }

        // send all non-zero entries columns
        for (long i=0; i<PI.b; i++) {
            if (curr_colsum[i] > 0) {
                bsp_send(t,
                         &i,
                         &curr_colsum[i],
                         sizeof(long)
                );
            }

            curr_colsum[i] = 0;
        }
    }

    bsp_sync();

    // determine number of messages received
    int nr_msgs;
    bsp_qsize(&nr_msgs, NULL);

    // for each message, increment the column
    for (long m=0; m<nr_msgs; m++) {
        bsp_size_t status;
        long col_local;
        long col_value;
        bsp_get_tag(&status, &col_local);
        bsp_move(&col_value, sizeof(long));
        colsum[col_local] += col_value;
    }
    
    free(start_indices);
    free(curr_indices);
    free(curr_colsum);

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

    for (bsp_pid_t t=0; t<=PI.p_last; t++) {
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