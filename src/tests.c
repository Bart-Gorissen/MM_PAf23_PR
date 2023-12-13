#include <stdlib.h>
#include <time.h>
#include <bsp.h>

#include "utility.h"
#include "CRSgraph.h"
#include "utility_parallel.h"
#include "pagerank_sequential.h"
#include "pagerank_parallel.h"
#include "tests.h"

/**
 * Checks whether the random function produces values from [0, N)
 * This is done by repeated sampling, until all numbers are seen
 * NOTE: numbers of values outside the range are only reported on termination
 * Return: iters < max_iters
*/
int random_test(long N, long max_iters, int V) {

    init_randomness(0);

    int* seen_randoms = (int*) calloc(N, sizeof(int));
    long nr_unique = 0;
    long nr_larger = 0;
    long iters = 0;

    while (nr_unique < N && iters < max_iters) {
        long val = urandom_rrange(N);
        if (val >= N) {
            nr_larger += 1;
        }
        else if (seen_randoms[val] == 0) {
            seen_randoms[val] = 1;
            nr_unique ++;
        }
        iters ++;
    }    

    if (V > 0 && iters < max_iters) {
        printf("Completed random range test in %ld iterations, found %ld numbers outside of range\n--------\n", iters, nr_larger);
    }
    else if (V > 0 && iters >= max_iters) {
        printf("Failed random range test in %ld iterations, found %ld numbers outside of range\n", iters, nr_larger);
        if (V > 3) {
            printf("[ ");
            for (long i=0; i<N; i++) {
                printf("%d ", seen_randoms[i]);
            }
            printf("]\n");
        }
        printf("--------\n");
    }

    free(seen_randoms);

    return iters < max_iters;
}

/**
 * Tests graph construction, sorting, column sum, and stochastic diagonal
 * Return: 1
*/
int graph_seq_test(long N, long L, double p, int V) {

    init_randomness(0);

    CRSGraph graph = CRSGraph_generate(N, N, L);
    long* colsum_pre = CRSGraph_colsum(graph);
    long* diagonal_pre = CRSGraph_stochastic_diagonal(graph);

    if (V > 10) {
        CRSGraph_print(graph);
    }
    if (V > 5) {
        printf("[ ");
        for (long i=0; i<graph.N; i++) {
            printf("%ld ", colsum_pre[i]);
        }
        printf("] Column Sum\n");
        printf("[ ");
        for (long i=0; i<graph.N; i++) {
            printf("%ld ", diagonal_pre[i]);
        }
        printf("] Stochastic Diagonal\n");
    }
    
    CRSGraph_sort(graph);
    long* colsum_post = CRSGraph_colsum(graph);
    long* diagonal_post = CRSGraph_stochastic_diagonal(graph);

    if (V > 10) {
        CRSGraph_print(graph);
    }
    if (V > 5) {
        printf("[ ");
        for (long i=0; i<graph.N; i++) {
            printf("%ld ", colsum_post[i]);
        }
        printf("] Column Sum\n");
        printf("[ ");
        for (long i=0; i<graph.N; i++) {
            printf("%ld ", diagonal_post[i]);
        }
        printf("] Stochastic Diagonal\n");
    }
    
    free(graph.rowsize);
    free(graph.colindex);
    free(colsum_pre);
    free(colsum_post);
    free(diagonal_pre);
    free(diagonal_post);

    if (V > 0) {
        printf("Completed Graph (sequential) test 0\n--------\n");
    }

    return 1;
}

/**
 * Prints the result of sorting CRS representation
 * Return: 1
*/
int graph_sort_seq_test(long N, long L, int V) {

    init_randomness(0);

    CRSGraph graph = CRSGraph_generate(N, N, L);
    long start = 0;

    if (V > 2) {
        start = 0;
        printf("[ ");
        for (long i=0; i<N; i++) {
            for (long j=0; j<graph.rowsize[i]; j++) {
                printf("%ld ", graph.colindex[start+j]);
            }
            if (i < N-1) printf("| ");
            start += graph.rowsize[i];
        }
        printf("]\n[ ");
        for (long i=0; i<N; i++) {
            printf("%ld ", graph.rowsize[i]);
        }
        printf("]\n");
    }
    
    CRSGraph_sort(graph);

    if (V > 2) {
        start = 0;
        printf("[ ");
        for (long i=0; i<N; i++) {
            for (long j=0; j<graph.rowsize[i]; j++) {
                printf("%ld ", graph.colindex[start+j]);
            }
            if (i < N-1) printf("| ");
            start += graph.rowsize[i];
        }
        printf("]\n[ ");
        for (long i=0; i<N; i++) {
            printf("%ld ", graph.rowsize[i]);
        }
        printf("]\n");
    }

    free(graph.rowsize);
    free(graph.colindex);

    if (V > 0) {
        printf("Completed Graph (sequential) test 1\n--------\n");
    }

    return 1;
}


/**
 * Tests pagerank_seq function and scaling of resulting vector.
*/
int pagerank_seq_test1(long N, long L, double p, double eps, int V) {

    init_randomness(0);

    CRSGraph graph = CRSGraph_generate(N, N, L);

    double* u_sln = pagerank_seq(graph, p, eps, 0, 0, 2);
    double u_len = 0;
    double u_len_post = 0;
    long nr_negative = 0;

    for (long i=0; i<graph.N; i++) {
        u_len += u_sln[i];
        if (u_sln[i] < 0) {
            nr_negative ++;
        }
    }
    for (long i=0; i<graph.N; i++) {
        u_len_post += u_sln[i] / u_len;
    }

    if (V > 5) {
        printf("[ ");
        for (long i=0; i<graph.N; i++) {
            printf("%.3f ", u_sln[i]);
        }
        printf("] Pagerank Solution\n");
        
        printf("[ ");
        for (long i=0; i<graph.N; i++) {
            printf("%.3f ", u_sln[i] / u_len);
        }
        printf("] Pagerank Solution (Scaled)\n");
    }
    if (V > 3) {
        printf("Length after scaling: %.3f\n", u_len_post);
        printf("Negative values: %ld\n", nr_negative);
    }

    free(graph.rowsize);
    free(graph.colindex);
    free(u_sln);

    if (V > 0) {
        printf("Completed Pagerank (sequential) test 1\n--------\n");
    }

    return 1;
}

/**
 * Compares running time and number of iterations between sorted and unsorted CRS
*/
int pagerank_seq_test2(long N, long L, double p, double eps, int V) {

    init_randomness(0);

    CRSGraph graph = CRSGraph_generate(N, N, L);
    clock_t start = clock();
    double* u_sln_pre = pagerank_seq(graph, p, eps, 0, 0, V);
    printf("Unsorted took time %f\n", ((double) clock()-start)/CLOCKS_PER_SEC);

    CRSGraph_sort(graph);
    start = clock();
    double* u_sln_post = pagerank_seq(graph, p, eps, 0, 0, V);
    printf("Sorted took time %f\n", ((double) clock()-start)/CLOCKS_PER_SEC);
    
    free(graph.rowsize);
    free(graph.colindex);
    free(u_sln_pre);
    free(u_sln_post);

    if (V > 0) {
        printf("Completed Pagerank (sequential) test 2\n--------\n");
    }

    return 1;
}

/**
 * Compares running time and number of iterations between standard and pre-inverted D
*/
int pagerank_seq_test3(long N, long L, double p, double eps, int V) {

    init_randomness(0);

    CRSGraph graph = CRSGraph_generate(N, N, L);
    clock_t start = clock();
    double* u_sln_pre = pagerank_seq(graph, p, eps, 0, 0, V);
    printf("Standard D took time %f\n", ((double) clock()-start)/CLOCKS_PER_SEC);

    start = clock();
    double* u_sln_post = pagerank_seq(graph, p, eps, 1, 0, V);
    printf("Inverted D took time %f\n", ((double) clock()-start)/CLOCKS_PER_SEC);
    
    free(graph.rowsize);
    free(graph.colindex);
    free(u_sln_pre);
    free(u_sln_post);

    if (V > 0) {
        printf("Completed Pagerank (sequential) test 3\n--------\n");
    }

    return 1;
}

/**
 * Compares iterations between fast and full residual computation
*/
int pagerank_seq_test4(long N, long L, double p, double eps, int V) {

    init_randomness(0);

    CRSGraph graph = CRSGraph_generate(N, N, L);
    clock_t start = clock();
    double* u_sln_pre = pagerank_seq(graph, p, eps, 0, 0, V);
    printf("Fast residual took time %f\n", ((double) clock()-start)/CLOCKS_PER_SEC);

    start = clock();
    double* u_sln_post = pagerank_seq(graph, p, eps, 0, 1, V);
    printf("Full residual took time %f\n", ((double) clock()-start)/CLOCKS_PER_SEC);
    
    free(graph.rowsize);
    free(graph.colindex);
    free(u_sln_pre);
    free(u_sln_post);

    if (V > 0) {
        printf("Completed Pagerank (sequential) test 4\n--------\n");
    }

    return 1;
}

/**
 * Tests pagerank_par_naive function.
*/
int pagerank_par_test1(long N, long L, double p, double eps, int V) {

    PARInfo PI = pagerank_par_init(N, 0);
    CRSGraph graph = CRSGraph_generate(N, PI.b_local, L);

    double* u_sln = pagerank_par(graph, p, eps, PI, 0, 0, 0, 0, V);
    double u_len = 0;
    double u_len_post = 0;
    long nr_negative = 0;

    for (long i=0; i<PI.b_local; i++) {
        u_len += u_sln[i];
        if (u_sln[i] < 0) {
            nr_negative ++;
        }
    }
    for (long i=0; i<PI.b_local; i++) {
        u_len_post += u_sln[i] / u_len;
    }

    if (V > 5) {
        printf("[ ");
        for (long i=0; i<PI.b_local; i++) {
            printf("%.3f ", u_sln[i]);
        }
        printf("] Pagerank Solution\n");
        
        printf("[ ");
        for (long i=0; i<PI.b_local; i++) {
            printf("%.3f ", u_sln[i] / u_len);
        }
        printf("] Pagerank Solution (Scaled)\n");
    }
    if (V > 3) {
        printf("Length after scaling: %.3f\n", u_len_post);
        printf("Negative values: %ld\n", nr_negative);
    }

    free(graph.rowsize);
    free(graph.colindex);
    free(u_sln);

    if (V > 0) {
        printf("Completed Pagerank (parallel naive) test 1\n--------\n");
    }

    return 1;
}

/**
 * Tests pagerank_par_gget function.
*/
int pagerank_par_test2(long N, long L, double p, double eps, int V) {

    PARInfo PI = pagerank_par_init(N, 0);
    CRSGraph graph = CRSGraph_generate(N, PI.b_local, L);

    double* u_sln = pagerank_par(graph, p, eps, PI, 0, 0, 0, 2, V);
    double u_len = 0;
    double u_len_post = 0;
    long nr_negative = 0;

    for (long i=0; i<PI.b_local; i++) {
        u_len += u_sln[i];
        if (u_sln[i] < 0) {
            nr_negative ++;
        }
    }
    for (long i=0; i<PI.b_local; i++) {
        u_len_post += u_sln[i] / u_len;
    }

    if (V > 5) {
        printf("[ ");
        for (long i=0; i<PI.b_local; i++) {
            printf("%.3f ", u_sln[i]);
        }
        printf("] Pagerank Solution\n");
        
        printf("[ ");
        for (long i=0; i<PI.b_local; i++) {
            printf("%.3f ", u_sln[i] / u_len);
        }
        printf("] Pagerank Solution (Scaled)\n");
    }
    if (V > 3) {
        printf("Length after scaling: %.3f\n", u_len_post);
        printf("Negative values: %ld\n", nr_negative);
    }

    free(graph.rowsize);
    free(graph.colindex);
    free(u_sln);

    if (V > 0) {
        printf("Completed Pagerank (parallel gget) test 2\n--------\n");
    }

    return 1;
}

/**
 * Tests different computations for D function.
*/
int pagerank_par_test3(long N, long L, double p, double eps, int V) {

    PARInfo PI = pagerank_par_init(N, 0);
    CRSGraph graph = CRSGraph_generate(N, PI.b_local, L);

    double* u_sln;
    u_sln = pagerank_par(graph, p, eps, PI, 0, 0, 0, 0, V);
    free(u_sln);
    u_sln = pagerank_par(graph, p, eps, PI, 0, 0, 1, 0, V);
    free(u_sln);
    u_sln = pagerank_par(graph, p, eps, PI, 0, 0, 2, 0, V);
    free(u_sln);
    u_sln = pagerank_par(graph, p, eps, PI, 0, 0, 3, 0, V);
    free(u_sln);

    free(graph.rowsize);
    free(graph.colindex);

    if (V > 0) {
        printf("Completed Pagerank (D computation) test 3\n--------\n");
    }

    return 1;
}