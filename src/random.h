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
 * @param weight number of nonzero coefficients
 */
void random_weighted_gf4_poly(gf4_poly_t * poly, size_t weight);

#endif //MDPC_GF4_RANDOM_H
