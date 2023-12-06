#include "utility.h"
#include "CRSgraph.h"
#include "utility_parallel.h"
#include "pagerank_parallel.h"

/**
 * Initializes random seed and parallel info
 * Return: PARInfo on BSP setup and data distribution
*/
PARInfo pagerank_par_init(long N, long seed) {
    PARInfo PI = setup_parinfo(N);
    init_randomness(seed + PI.s);
    return PI;
}

/**
 * Parallel pagerank implementation
 * NOTE: naive: communicates full x to all processors during pGD^{-1}x
 * Return: local part of pagerank solution
*/
double* pagerank_par_naive(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V) {
    // setup communication things
    double* p_vec_par = (double*) malloc(PI.P * sizeof(double)); // vector for BSP use with P entries
    double* r_vec_par = (double*) malloc(PI.b_local * sizeof(double)); // vector for BSP use with b_local entries
    bsp_push_reg(p_vec_par, PI.P * sizeof(double));
    bsp_push_reg(r_vec_par, PI.b_local * sizeof(double));
    bsp_sync();

    // compute the matrix D
    long* D_diag = parallel_colsum_bigpull(graph, PI);
    change_01(D_diag, PI.b_local); // change 0 entries to 1 entries

    // make stochastic vector u
    double* u_vec = (double*) malloc(PI.b_local * sizeof(double));
    double local_sum = generate_urandom01(u_vec, PI.b_local);
    double global_1norm = parallel_sum(local_sum, p_vec_par, PI);

    // compute initial residual r = e - (I-pGD^{-1})u = e - u + pGD^{-1} u    
    // first scale all u entries w.r.t D (r = D^{-1} u)
    for (long i=0; i<PI.b_local; i++) {
        r_vec_par[i] = u_vec[i] / D_diag[i];
    }

    // compute r = p G_s (D^{-1} u) [global operation]
    //bsp_sync();
    parallel_pGy_bigpull(graph, p, r_vec_par, r_vec_par, PI);

    // compute r = e - u + (pG_s (D^{-1} u))
    for (long i=0; i<PI.b_local; i++) {
        r_vec_par[i] = 1 - u_vec[i] + r_vec_par[i];
    }

    bsp_sync();

    double residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
    long iters = 0;

    // pagerank loop 
    while (residual_norm >= eps) {
        // compute next u = u + r
        // compute next r = pGD^{-1} r

        // first set u = u + r and scale all r entries w.r.t D: r = D^{-1} r
        for (long i=0; i<PI.b_local; i++) {
            u_vec[i] += r_vec_par[i];
            r_vec_par[i] /= D_diag[i];
        }

        // compute r = p G_s . (D^{-1} r)
        //bsp_sync();
        parallel_pGy_bigpull(graph, p, r_vec_par, r_vec_par, PI);
        residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
        iters ++;

        if (V > 4) printf("Current norm: %f. Expected next norm: %f\n", residual_norm, p*residual_norm);   
    }

    if (V > 2) printf("Took %ld iterations\n", iters);

    bsp_pop_reg(p_vec_par);
    bsp_pop_reg(r_vec_par);

    free(D_diag);
    free(p_vec_par);
    free(r_vec_par);

    return u_vec;
}