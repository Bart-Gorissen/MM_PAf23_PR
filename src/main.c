#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <bsp.h>

#include "utility.h"
#include "tests.h"

PARAMS PRMS;

int seq_tests() {
    random_test(PRMS.N, PRMS.N * PRMS.N, PRMS.V);
    graph_seq_test(PRMS.N, PRMS.L, PRMS.p, PRMS.V);
    graph_sort_seq_test(PRMS.N, PRMS.L, PRMS.V);
    pagerank_seq_test1(PRMS.N, PRMS.L, PRMS.p, PRMS.eps, PRMS.V);
    pagerank_seq_test2(PRMS.N, PRMS.L, PRMS.p, PRMS.eps, PRMS.V);
    pagerank_seq_test3(PRMS.N, PRMS.L, PRMS.p, PRMS.eps, PRMS.V);

    exit(EXIT_SUCCESS);
}

void par_test_outer() {
    bsp_begin(PRMS.P);

    pagerank_par_test1(PRMS.N, PRMS.L, PRMS.p, PRMS.eps, PRMS.V);

    bsp_end();
}




int main(int argc, char **argv) {
    bsp_init(par_test_outer, argc, argv);

    int has_parsed = parse_input(&PRMS, argc, argv);
    if (has_parsed == 0) {
        exit(EXIT_FAILURE);
    }

    seq_tests();

    //par_test_outer();

    exit(EXIT_SUCCESS);
}