#include <stdio.h>
#include <bsp.h>

#include "utility.h"
#include "CRSgraph.h"
#include "pagerank_sequential.h"

/**
 * Computes pagerank solution
 * pre_invert_D {0,1} : computes D^{-1} in memory
 * u_choice {0,1,...,N}: (0 : random), (i : e_i)
 * full_res {0,1} : ( 0 : computes r update from r ), ( 1 : computes r from u )
 * Verbosity level: (0 - nothing) (1 - number of iterations) (2 - norm residual)
 * Return: solution u of system (I - pGD^{-1})u = e with residual norm less than eps
*/
double* pagerank_seq(CRSGraph graph, double p, double eps, int pre_invert_D, int u_choice, int full_res, int V) {

    double time_PR_total = bsp_time();
    double time_D_comm = bsp_time();
    double time_pGDx_comm;

    // compute the matrix D and prepare vectors u and r
    long* D_diag = CRSGraph_stochastic_diagonal(graph);
    double* D_diag_inv = NULL;

    // pre-invert if needed
    if (pre_invert_D != 0) {
        D_diag_inv = inverse(D_diag, graph.N);
    }

    time_D_comm = bsp_time() - time_D_comm;

    double* u_vec;

    if (u_choice == 0) {
        u_vec = generate_stochastic(graph.N);
    }
    else if (u_choice < 0) {
        u_vec = generate_vector_filled(1.0 / graph.N, graph.N);
    }
    else {
        u_vec = generate_ei(u_choice-1, graph.N);
    }

    double* r_vec = (double*) malloc(graph.N * sizeof(double));
    double* r_aux = (double*) malloc(graph.N * sizeof(double));

    time_pGDx_comm = bsp_time();
    
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

    time_pGDx_comm = bsp_time() - time_pGDx_comm;

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

    time_PR_total = bsp_time() - time_PR_total;

    // write output to file
    if (V == -8) {
        char filename[42];
        sprintf(filename, "out/output_%ld.out", graph.N);
        FILE* file = fopen(filename, "a+");
        fprintf(file, "%ld, %lf, %lf, %lf\n", iters, time_D_comm, time_pGDx_comm, time_PR_total);
        fclose(file);
    }

    return u_vec;
}