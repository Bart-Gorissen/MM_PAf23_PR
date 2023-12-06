#ifndef UTILITY_H
#define UTILITY_H

#include <stdlib.h>
#include <math.h>
#include <bsp.h>

// Parameter struct
typedef struct {
    long N;         // number of nodes
    long L;         // inlink number
    double p;       // follow probability
    long P;    // number of processors
    double eps;     // precision
    int V;          // verbosity
    long seed; // starting seed
} PARAMS;

/**
 * Reads input from arguments into PRMS
 * Return: (0 - failure), (1 - success)
*/
int parse_input(PARAMS* PRMS, int argc, char **argv);

/**
 * Prepares pseudo-random generator (sets seed)
 * To be called once at initialization of program
 * Return: void
*/
void init_randomness(long seed);

/**
 * Returns a random value from [0, right_bound)
 * Expects right_bound > 0
 * Return: [0,right_bound)
*/
long urandom_rrange(long right_bound);

/**
 * Returns a random value from [0,1]
 * Return: [0,1]
*/
double urandom_01();

/**
 * Sorts array arr between indices start and end (inclusive)
 * Return: void
*/
void quicksort(long* arr, long start, long end);

/**
 * Fills vector result with N entries in [0,1]
 * Return: sum of entries
*/
double generate_urandom01(double* result, long N);

/**
 * Generates stochastic vector of length N
 * Return: stochastic vector
*/
double* generate_stochastic(long N);

/**
 * Changes 0 entries in arr of length N to 1
 * Return: void
*/
void change_01(long* arr, long N);

/**
 * Computes the sum of squares of N entries of u_vec
 * Return: sum of squares
*/
double sq_sum(double* u_vec, long N);

/**
 * Computes the (Euclidean) norm of N entries of u_vec
 * Return: Euclidean norm u_vec
*/
double norm(double* u_vec, long N);

/**
 * Computes the inverse of full-rank (integer) diagonal matrix D
 * Return: D^{-1}
*/
double* inverse(long* D, long N);

#endif /* UTILITY_H */