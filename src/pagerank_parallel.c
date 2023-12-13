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
 * pre_invert_D {0,1}: pre-computing the inverse of D
 * u_choice {0,1,...,N}: (0 : random), (i : e_i)
 * D_comm_choice {0,1,2}: (0 : full-broadcast), (1 : P-round broadcast), (2 : message queue)
 * pGx_comm_choice {0,1,2}: (0 : full-broadcast of u), (1 : P-round broadcast), (2 : quick-get), (3 : mapped-get)
 * Return: local part of pagerank solution
*/
double* pagerank_par(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int u_choice, int D_comm_choice, int pGx_comm_choice, int V) {
    // setup communication vectors (whether these are used depends on the parameters)
    double* p_vec_par = (double*) malloc(PI.P * sizeof(double)); // vector for BSP use with P entries
    double* r_vec_par = (double*) malloc(PI.b_local * sizeof(double)); // vector for BSP use with b_local entries
    long* D_vec_par;
    bsp_size_t tagsize;
    switch(D_comm_choice) {
        case 0:
            D_vec_par = (long*) calloc(graph.N, sizeof(long));
            bsp_push_reg(D_vec_par, graph.N * sizeof(long));
            break;

        case 1:
            D_vec_par = (long*) calloc(PI.b, sizeof(long));
            bsp_push_reg(D_vec_par, PI.b * sizeof(long));
            break;
        
        case 2:
            D_vec_par = (long*) malloc(0 * sizeof(long));
            bsp_push_reg(D_vec_par, 0 * sizeof(long));
            tagsize = 0;
            bsp_set_tagsize(&tagsize);
            break;

        case 3:
            D_vec_par = (long*) malloc(0 * sizeof(long));
            bsp_push_reg(D_vec_par, 0 * sizeof(long));
            tagsize = sizeof(long);
            bsp_set_tagsize(&tagsize);
            break;
        
        default:
            break;
    }
    bsp_push_reg(p_vec_par, PI.P * sizeof(double));
    bsp_push_reg(r_vec_par, PI.b_local * sizeof(double));

    bsp_sync();

    // generate vector u
    double* u_vec;

    if (u_choice == 0) {
        // random vector
        u_vec = (double*) malloc(PI.b_local * sizeof(double));
        double local_sum = generate_urandom01(u_vec, PI.b_local);
        double global_1norm = parallel_sum(local_sum, p_vec_par, PI);
        
        // normalize u
        for (long i=0; i<PI.b_local; i++) {
            u_vec[i] /= global_1norm;
        }
    }
    else {
        // vector e_{u_choice}
        u_vec = (double*) calloc(PI.b_local, sizeof(double));

        // check if 1-entry is on this processor
        if ( (u_choice-1) / PI.b == PI.s ) {
            u_vec[(u_choice-1) % PI.b] = 1;
        }
    }

    // compute the matrix D
    long* D_diag = parallel_colsum(graph, PI, D_vec_par, D_comm_choice);    

    change_01(D_diag, PI.b_local); // change 0 entries to 1 entries
    double* D_diag_inv = NULL;

    // pre-invert if needed
    if (pre_invert_D != 1) {
        D_diag_inv = inverse(D_diag, PI.b_local);
    }

    // compute initial residual r = e - (I-pGD^{-1})u = e - u + pGD^{-1} u    
    // first scale all u entries w.r.t D (r = D^{-1} u)
    if (pre_invert_D == 0) {
        for (long i=0; i<PI.b_local; i++) {
            r_vec_par[i] = u_vec[i] / D_diag[i];
        }
    }
    else {
        for (long i=0; i<PI.b_local; i++) {
            r_vec_par[i] = u_vec[i] * D_diag_inv[i];
        }
    }
    
    // compute r = p G_s (D^{-1} u) [global operation]
    parallel_pGx(graph, p, r_vec_par, r_vec_par, PI, pGx_comm_choice);

    // compute r = e - u + (pG_s (D^{-1} u))
    for (long i=0; i<PI.b_local; i++) {
        r_vec_par[i] = 1 - u_vec[i] + r_vec_par[i];
    }

    double residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
    long iters = 0;

    if (V > 2) printf("Expecting %ld iterations\n", (long) (log(eps / residual_norm) / log(p)));

    // pagerank loop 
    while (residual_norm >= eps) {
        // compute next u = u + r
        // compute next r = pGD^{-1} r

        // first set u = u + r and scale all r entries w.r.t D: r = D^{-1} r
        if (pre_invert_D == 0) {
            for (long i=0; i<PI.b_local; i++) {
                u_vec[i] += r_vec_par[i];
                r_vec_par[i] /= D_diag[i];
            }
        }
        else {
            for (long i=0; i<PI.b_local; i++) {
                u_vec[i] += r_vec_par[i];
                r_vec_par[i] *= D_diag_inv[i];
            }
        }

        // compute r = p G_s . (D^{-1} r)
        parallel_pGx(graph, p, r_vec_par, r_vec_par, PI, pGx_comm_choice);
        
        residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
        iters ++;

        if (V > 4) printf("Current norm: %f. Expected next norm: %f\n", residual_norm, p*residual_norm);   
    }

    if (V > 2) printf("Took %ld iterations\n", iters);
    
    bsp_pop_reg(p_vec_par);
    bsp_pop_reg(r_vec_par);
    bsp_pop_reg(D_vec_par);
    free(p_vec_par);
    free(r_vec_par);
    free(D_vec_par);
    free(D_diag);
    free(D_diag_inv);

    return u_vec;
}

// /**
//  * [Depricated] Parallel pagerank implementation
//  * NOTE: naive: communicates full x to all processors during pGD^{-1}x
//  * Return: local part of pagerank solution
// */
// double* pagerank_par_naive(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V) {
//     // setup communication things
//     double* p_vec_par = (double*) malloc(PI.P * sizeof(double)); // vector for BSP use with P entries
//     double* r_vec_par = (double*) malloc(PI.b_local * sizeof(double)); // vector for BSP use with b_local entries
//     bsp_push_reg(p_vec_par, PI.P * sizeof(double));
//     bsp_push_reg(r_vec_par, PI.b_local * sizeof(double));
//     bsp_sync();

//     // compute the matrix D
//     long* D_diag = parallel_colsum_bigpull(graph, PI);
//     change_01(D_diag, PI.b_local); // change 0 entries to 1 entries

//     // make stochastic vector u
//     double* u_vec = (double*) malloc(PI.b_local * sizeof(double));
//     double local_sum = generate_urandom01(u_vec, PI.b_local);
//     double global_1norm = parallel_sum(local_sum, p_vec_par, PI);

//     // compute initial residual r = e - (I-pGD^{-1})u = e - u + pGD^{-1} u    
//     // first scale all u entries w.r.t D (r = D^{-1} u)
//     for (long i=0; i<PI.b_local; i++) {
//         r_vec_par[i] = u_vec[i] / D_diag[i];
//     }

//     // compute r = p G_s (D^{-1} u) [global operation]
//     //bsp_sync();
//     parallel_pGy_bigpull(graph, p, r_vec_par, r_vec_par, PI);

//     // compute r = e - u + (pG_s (D^{-1} u))
//     for (long i=0; i<PI.b_local; i++) {
//         r_vec_par[i] = 1 - u_vec[i] + r_vec_par[i];
//     }

//     bsp_sync();

//     double residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
//     long iters = 0;

//     // pagerank loop 
//     while (residual_norm >= eps) {
//         // compute next u = u + r
//         // compute next r = pGD^{-1} r

//         // first set u = u + r and scale all r entries w.r.t D: r = D^{-1} r
//         for (long i=0; i<PI.b_local; i++) {
//             u_vec[i] += r_vec_par[i];
//             r_vec_par[i] /= D_diag[i];
//         }

//         // compute r = p G_s . (D^{-1} r)
//         //bsp_sync();
//         parallel_pGy_bigpull(graph, p, r_vec_par, r_vec_par, PI);
//         residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
//         iters ++;

//         if (V > 4) printf("Current norm: %f. Expected next norm: %f\n", residual_norm, p*residual_norm);   
//     }

//     if (V > 2) printf("Took %ld iterations\n", iters);

//     bsp_pop_reg(p_vec_par);
//     bsp_pop_reg(r_vec_par);

//     free(D_diag);
//     free(p_vec_par);
//     free(r_vec_par);

//     return u_vec;
// }

// /**
//  * [Depricated] Parallel pagerank implementation
//  * NOTE: gget: communicates only elements needed to processors during pGD^{-1}x, but possibly with multiplicity
//  * Return: local part of pagerank solution
// */
// double* pagerank_par_gget(CRSGraph graph, double p, double eps, PARInfo PI, int pre_invert_D, int V) {
//     // setup communication things
//     double* p_vec_par = (double*) malloc(PI.P * sizeof(double)); // vector for BSP use with P entries
//     double* r_vec_par = (double*) malloc(PI.b_local * sizeof(double)); // vector for BSP use with b_local entries
//     bsp_push_reg(p_vec_par, PI.P * sizeof(double));
//     bsp_push_reg(r_vec_par, PI.b_local * sizeof(double));
//     bsp_sync();

//     // compute the matrix D
//     long* D_diag = parallel_colsum_bigpull(graph, PI);
//     change_01(D_diag, PI.b_local); // change 0 entries to 1 entries

//     // make stochastic vector u
//     double* u_vec = (double*) malloc(PI.b_local * sizeof(double));
//     double local_sum = generate_urandom01(u_vec, PI.b_local);
//     double global_1norm = parallel_sum(local_sum, p_vec_par, PI);

//     // compute initial residual r = e - (I-pGD^{-1})u = e - u + pGD^{-1} u    
//     // first scale all u entries w.r.t D (r = D^{-1} u)
//     for (long i=0; i<PI.b_local; i++) {
//         r_vec_par[i] = u_vec[i] / D_diag[i];
//     }

//     // compute r = p G_s (D^{-1} u) [global operation]
//     //bsp_sync();
//     parallel_pGy_gget(graph, p, r_vec_par, r_vec_par, PI);

//     // compute r = e - u + (pG_s (D^{-1} u))
//     for (long i=0; i<PI.b_local; i++) {
//         r_vec_par[i] = 1 - u_vec[i] + r_vec_par[i];
//     }

//     bsp_sync();

//     double residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
//     long iters = 0;

//     // pagerank loop 
//     while (residual_norm >= eps) {
//         // compute next u = u + r
//         // compute next r = pGD^{-1} r

//         // first set u = u + r and scale all r entries w.r.t D: r = D^{-1} r
//         for (long i=0; i<PI.b_local; i++) {
//             u_vec[i] += r_vec_par[i];
//             r_vec_par[i] /= D_diag[i];
//         }

//         // compute r = p G_s . (D^{-1} r)
//         //bsp_sync();
//         parallel_pGy_gget(graph, p, r_vec_par, r_vec_par, PI);
//         residual_norm = parallel_norm(r_vec_par, p_vec_par, PI);
//         iters ++;

//         if (V > 4) printf("Current norm: %f. Expected next norm: %f\n", residual_norm, p*residual_norm);   
//     }

//     if (V > 2) printf("Took %ld iterations\n", iters);

//     bsp_pop_reg(p_vec_par);
//     bsp_pop_reg(r_vec_par);

//     free(D_diag);
//     free(p_vec_par);
//     free(r_vec_par);

//     return u_vec;
// }