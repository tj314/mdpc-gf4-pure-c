#ifndef MDPC_GF4_RANDOM_H
#define MDPC_GF4_RANDOM_H

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "gf4_poly.h"
#include "gf4.h"


/**
 * Initialize random using srand with current time.
 *
 * Ensures that srand is called at most once during each program execution. This function is not thread safe.
 * This function doesn't need to be called explicitly.
 */
void random_init();

/**
 * Force another call to srand using current time.
 */
void random_force_reseed();

/**
 * Return an unsigned integer within given range.
 *
 * May also call random_init().
 * This function uses rand() to generate random values.
 *
 * @param low_bound_inclusive inclusive low bound
 * @param top_bound_inclusive inclusive top bound
 * @return an unsigned integer i s. t. low_bound_inclusive <= i <= top_bound_inclusive
 */
size_t random_from_range(size_t low_bound_inclusive, size_t top_bound_inclusive);

/**
 * Generate a polynomial of given weight.
 *
 * May also call random_init().
 *
 * @param poly output polynomial, it must be allocated in advance
 * @param size maximum size of the polynomial, size <= polynomial->capacity
 * @param weight number of nonzero coefficients, weight <= size
 */
void random_weighted_gf4_poly(gf4_poly_t * poly, size_t size, size_t weight);


/**
 * Generate a polynomial of given weight such that at least weight/2 ones are placed exactly distance apart.
 *
 * @param poly
 * @param size
 * @param weight
 * @param distance
 */
void random_weighted_gf4_poly_pairs_of_ones(gf4_poly_t * poly, size_t size, size_t weight, size_t distance);

void random_weighted_gf4_poly_pairs_of_one_alpha(gf4_poly_t * poly, size_t size, size_t weight, size_t distance);

#endif //MDPC_GF4_RANDOM_H
