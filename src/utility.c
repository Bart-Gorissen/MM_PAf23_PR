#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <bsp.h>

#include "utility.h"

/**
 * Reads input from arguments into PRMS
 * Return: (0 - failure), (1 - success)
*/
int parse_input(PARAMS* PRMS, int argc, char **argv) {

    PRMS->N=1;       // number of nodes in G
    PRMS->L=1;       // inlink number
    PRMS->p=0.5;     // chance follow
    PRMS->P=1;       // number of processors
    PRMS->eps = 0.000001; // 10^{-6}
    PRMS->V=0;       // verbosity level
    PRMS->seed = 42069314; // starting seed

    int input_type = 0; // 0: interactive, 1: arguments, 2: arguments(flag), 3: file
    char* input_file = NULL;

    // check which input type we are using
    if (argc == 1) {
        input_type = 0;
    }
    else { // argc > 1
        long split_index = -1;
        for (long j=0; j<strlen(argv[1]); j++) {
            if (argv[1][j] == '=') {
                split_index = j;
                break;
            }
        }

        if (split_index < 0) { // no file "flag"
            input_type = 1;
        }
        else { // possible file "flag"
            char* head = (char*) calloc(split_index + 1, sizeof(char));
            char* tail = (char*) calloc(strlen(argv[1]) - split_index, sizeof(char));

            for (long j=0; j<split_index; j++) {
                head[j] = argv[1][j];
            }
            for (long j=0; j<strlen(argv[1]) - split_index - 1; j++) {
                tail[j] = argv[1][j + split_index + 1];
            }
            
            if (strcmp("-f", head) == 0) {
                input_type = 3;
                input_file = tail;
            }
            else {
                input_type = 2;
                free(tail);
            }
            free(head);
        }
    }

    __uint8_t is_valid = 1;

    switch (input_type)
    {
    case 0: // interactive
        printf("Number of nodes in G: ");
        fflush(stdout);
        is_valid = is_valid && scanf("%ld",&PRMS->N);
        printf("Inlink number: ");
        fflush(stdout);
        is_valid = is_valid && scanf("%ld",&PRMS->L);
        printf("Follow probability: ");
        fflush(stdout);
        is_valid = is_valid && scanf("%lf",&PRMS->p);
        printf("Number of processors: ");
        fflush(stdout);
        is_valid = is_valid && scanf("%ld",&PRMS->P);
        printf("Verbosity: ");
        fflush(stdout);
        is_valid = is_valid && scanf("%d",&PRMS->V);
        break;
    
    case 1: // arguments
        if (argc > 1) {
            PRMS->N = atol(argv[1]);
        }
        if (argc > 2) {
            PRMS->L = atol(argv[2]);
        }
        if (argc > 3) {
            PRMS->p = atof(argv[3]);
        }
        if (argc > 4) {
            PRMS->P = atol(argv[4]);
        }
        if (argc > 5) {
            PRMS->V = atol(argv[5]);
        }
        if (argc > 6) {
            printf("Ignoring sixth parameter and beyond\n");
        }
        break;
    
    case 2: // arguments (flag)
        for (int i=1; i<argc; i++) {
            long split_index = -1;
            for (long j=0; j<strlen(argv[i]); j++) {
                if (argv[i][j] == '=') {
                    split_index = j;
                    break;
                }
            }

            if (split_index < 0) {
                printf("No flag found on argument %d: %s\n", i-1, argv[i]);
                printf("Usage: pagerank_seq <N> <L> <p> <P> <V>\n");
                return 0;
            }

            // if split_index >= 0, then we are dealing with flag arguments
            char* head = (char*) calloc(split_index + 1, sizeof(char));
            char* tail = (char*) calloc(strlen(argv[i]) - split_index, sizeof(char));

            for (long j=0; j<split_index; j++) {
                head[j] = argv[i][j];
            }
            for (long j=0; j<strlen(argv[i]) - split_index - 1; j++) {
                tail[j] = argv[i][j + split_index + 1];
            }

            if (strcmp("-N", head) == 0) {
                PRMS->N = atol(tail);
            }
            else if (strcmp("-L", head) == 0) {
                PRMS->L = atol(tail);
            }
            else if (strcmp("-p", head) == 0) {
                PRMS->p = atof(tail);
            }
            else if (strcmp("-P", head) == 0) {
                PRMS->P = atol(tail);
            }
            else if (strcmp("-V", head) == 0) {
                PRMS->V = atol(tail);
            }
            else {
                is_valid = 0;
            }

            free(head);
            free(tail);
        }
        break;

    case 3: // case
        FILE* file = fopen(input_file, "r");
        if (file == NULL) {
            printf("Invalid file name %s\n", input_file);
            printf("Usage: pagerank_seq <N> <L> <p> <P> <V>\n");
            free(input_file);
            return 0;
        }

        int read_args = fscanf(file, "%ld\n%ld\n%lf\n%ld\n%d", &PRMS->N, &PRMS->L, &PRMS->p, &PRMS->P, &PRMS->V);
        if (PRMS->V > 0) {
            switch (read_args)
            {
            case 0:
                printf("Using defaults for N, L, p, P, V\n");
                break;
            case 1:
                printf("Using defaults for L, p, P, V\n");
                break;
            case 2:
                printf("Using defaults for p, P, V\n");
                break;
            case 3:
                printf("Using defaults for P, V\n");
                break;
            case 4:
                printf("Using defaults for V\n");
                break;
            case 5:
                printf("All param sixth parameter and beyond\n");
            
            default:
                printf("More than 5 parameters read\n");
                break;
            }
        }

        fclose(file);
        free(input_file);

        break;
    
    default: //case
        printf("Error: invalid input type\n");
        printf("Usage: pagerank_seq <N> <L> <p> <P> <V>\n");
        return 0;
        break;
    }

    // check validity of inputc
    if (is_valid) {
        if (PRMS->N < 0) {
            printf("Error: specify non-negative N\n");
            is_valid = 0;
        }
        if (PRMS->L < 0) {
            printf("Error: specify non-negative L\n");
            is_valid = 0;
        }
        if (PRMS->p < 0 || PRMS->p > 1) {
            printf("Error: specify probability p in range [0,1]\n");
            is_valid = 0;
        }
        if (PRMS->P < 1 || PRMS->P > bsp_nprocs()) {
            printf("Error: specify processors P in range [1,%d]\n", bsp_nprocs());
            is_valid = 0;
        }
    }

    if (is_valid == 0) {
        printf("Usage: pagerank_seq <N> <L> <p> <P> <V>\n");
        return 0;
    }
    
    if (PRMS->V > 0) {
        printf("Using parameters: [ N=%ld, L=%ld, p=%lf, P=%ld, V=%d ]\n--------\n", PRMS->N, PRMS->L, PRMS->p, PRMS->P, PRMS->V);
    }

    return 1;
}

/**
 * Prepares pseudo-random generator
 * To be called once at initialization of program
 * Return: void
*/
void init_randomness(long seed) {
    //srandom(seed);
    srand(seed);
}

/**
 * Returns a random value from [0, right_bound)
 * Expects right_bound > 0
 * Return: [0,right_bound)
*/
long urandom_rrange(long right_bound) {
    long res, div = RAND_MAX / (right_bound);

    do {
        //res = random() / div;
        res = rand() / div;
    } while (res >= right_bound);

    return res;
    
    //this is sometimes a bit over the right_bound
    //return random() / (RAND_MAX / right_bound);
}

/**
 * Returns a random value from [0,1]
 * Return: [0,1]
*/
double urandom_01() {
    //return (double) random() / (double) RAND_MAX;
    return (double) rand() / (double) RAND_MAX;
}

/**
 * Quicksort implementation
 * Return: void
*/
void quicksort(long* arr, long start, long end) {
    if (start >= end) {
        return;
    }

    long pivot = arr[end];
    long i = start-1;

    for (int j=start; j<end; j++) {
        if (arr[j] <= pivot) {
            i ++;
            long aux = arr[j];
            arr[j] = arr[i];
            arr[i] = aux;
        }
    }

    i ++;

    long aux = arr[end];
    arr[end] = arr[i];
    arr[i] = aux;

    quicksort(arr, start, i-1);
    quicksort(arr, i+1, end);
}

/**
 * Fills vector result with N entries in [0,1]
 * Return: sum of entries
*/
double generate_urandom01(double* result, long N) {
    double sum = 0;

    for (long i=0; i<N; i++) {
        result[i] = urandom_01();
        sum += result[i];
    }

    return sum;
}

/**
 * Generates stochastic vector of length N
 * Return: stochastic vector
*/
double* generate_stochastic(long N) {
    double* result = (double*) malloc(N * sizeof(double));
    double sum = generate_urandom01(result, N);
    
    for (long i=0; i<N; i++) {
        result[i] /= sum;
    }

    return result;
}

/**
 * Changes 0 entries in arr of length N to 1
 * Return: void
*/
void change_01(long* arr, long N) {
    for (long i=0; i<N; i++) {
        if (arr[i] == 0) {
            arr[i] ++;
        }
    }
}

/**
 * Computes the sum of squares of N entries of u_vec
 * Return: sum of squares
*/
double sq_sum(double* u_vec, long N) {
    double sq_sum = 0;
    for (long i=0; i<N; i++) {
        sq_sum += u_vec[i] * u_vec[i];
    }
    return sq_sum;
}

/**
 * Computes the (Euclidean) norm of N entries of u_vec
 * Return: Euclidean norm u_vec
*/
double norm(double* u_vec, long N) {
    return sqrt(sq_sum(u_vec, N));
}

/**
 * Computes the inverse of full-rank (integer) diagonal matrix D
 * Return: D^{-1}
*/
double* inverse(long* D, long N) {
    double* result = malloc(N * sizeof(double));

    for (long i=0; i<N; i++) {
        result[i] = (double) 1 / (double) D[i];
    }

    return result;
}

/**
 * Generates vector of length N with all entries a
 * Return: a\vec{1}
*/
double* generate_vector_filled(double a, long N) {
    double* result = malloc(N * sizeof(double));
    for (long i=0; i<N; i++) {
        result[i] = a;
    }
    return result;
}

/**
 * Generates (stochastic) e_i vector
 * Return: e_i
*/
double* generate_ei(long i, long N) {
    double* result = calloc(N, sizeof(double));
    result[i] = 1;
    return result;
}