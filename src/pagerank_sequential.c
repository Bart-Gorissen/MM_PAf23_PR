#include <stdio.h>

#include "utility.h"
#include "CRSgraph.h"
#include "pagerank_sequential.h"

/**
 * Computes pagerank solution
 * pre_invert_D {0,1} : computes D^{-1} in memory
 * full_res {0,1} : ( 0 : computes r update from r ), ( 1 : computes r from u )
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq(CRSGraph graph, double p, double eps, int pre_invert_D, int full_res, int V) {

    // compute the matrix D and prepare vectors u and r
    long* D_diag = CRSGraph_stochastic_diagonal(graph);
    double* D_diag_inv = NULL;

    // pre-invert if needed
    if (pre_invert_D != 0) {
        D_diag_inv = inverse(D_diag, graph.N);
    }

    double* u_vec = generate_stochastic(graph.N);
    double* r_vec = (double*) malloc(graph.N * sizeof(double));
    double* r_aux = (double*) malloc(graph.N * sizeof(double));
    
    // start by computing the initial residual r = e - (I-pGD^{-1})u = e - u + pGD^{-1} u
    // first scale all u entries w.r.t D
    if (pre_invert_D == 0) {
        for (long i=0; i<graph.N; i++) {
            r_aux[i] = u_vec[i] / D_diag[i];
        }
    }
    else {
        for (long i=0; i<graph.N; i++) {
            r_aux[i] = u_vec[i] * D_diag_inv[i];
        }
    }

    long start = 0;

    for (long i=0; i<graph.N; i++) {

        // compute GD^{-1} u
        r_vec[i] = 0;

        for (long j=0; j<graph.rowsize[i]; j++) {
            r_vec[i] += r_aux[graph.colindex[start + j]];
        }

        // compute e - u + pGD^{-1} u per row
        r_vec[i] = 1 - u_vec[i] + (p * r_vec[i]);
        start += graph.rowsize[i];
    }

    long iters = 0;
    double residual_norm = norm(r_vec, graph.N);

    if (V > 2) printf("Expecting %ld iterations\n", (long) (log(eps / residual_norm) / log(p)));

    // pagerank loop 
    while (residual_norm >= eps) {

        // compute next u = u + r
        // compute next r_aux = pGD^{-1} r

        // first set u = u + r and scale all r entries w.r.t D
        if (full_res == 0) {
            if (pre_invert_D == 0) {
                for (long i=0; i<graph.N; i++) {
                    u_vec[i] += r_vec[i];
                    r_vec[i] /= D_diag[i];
                }
            }
            else {
                for (long i=0; i<graph.N; i++) {
                    u_vec[i] += r_vec[i];
                    r_vec[i] *= D_diag_inv[i];
                }
            }

            long start = 0;

            for (long i=0; i<graph.N; i++) {
                r_aux[i] = 0;

                for (long j=0; j<graph.rowsize[i]; j++) {
                    r_aux[i] += r_vec[graph.colindex[start + j]];
                }    

                r_aux[i] *= p;
                start += graph.rowsize[i];
            }

            // update values by swapping auxiliary and main pointers for r
            double* swap_tmp = r_vec;
            r_vec = r_aux;
            r_aux = swap_tmp;
        }
        else {
            // first set u = u + r and scale all new u entries w.r.t D
            if (pre_invert_D == 0) {
                for (long i=0; i<graph.N; i++) {
                    u_vec[i] += r_vec[i];
                    r_aux[i] = u_vec[i] / D_diag[i];
                }
            }
            else {
                for (long i=0; i<graph.N; i++) {
                    u_vec[i] += r_vec[i];
                    r_aux[i] = u_vec[i] * D_diag_inv[i];
                }
            }
            
            long start = 0;

            for (long i=0; i<graph.N; i++) {

                // compute GD^{-1} u
                r_vec[i] = 0;

                for (long j=0; j<graph.rowsize[i]; j++) {
                    r_vec[i] += r_aux[graph.colindex[start + j]];
                }

                // compute e - u + pGD^{-1} u per row
                r_vec[i] = 1 - u_vec[i] + (p * r_vec[i]);
                start += graph.rowsize[i];
            }
        }

        residual_norm = norm(r_vec, graph.N);

        if (V > 4) printf("Current norm: %f. Expected next norm: %f\n", residual_norm, p*residual_norm);
        iters ++;
    }

    if (V > 2) printf("Took %ld iterations\n", iters);

    free(D_diag);
    free(D_diag_inv);
    free(r_vec);
    free(r_aux);

    return u_vec;
}