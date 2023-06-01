/**
 *  @file   randoms.h
 *  @brief  Random functions.
 *  @author Tomáš Vavro
 *  @date   2023-05-12
 ***********************************************/

/*
 This file is part of QC-MDPC McEliece over GF(4) implementation.
 Copyright (C) 2023 Tomáš Vavro

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MDPC_GF4_RANDOM_H
#define MDPC_GF4_RANDOM_H

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "gf4_poly.h"
#include "gf4.h"


/**
 * @brief Initialize random using srand with current time.
 *
 * Ensures that srand is called at most once during each program execution. This function is not thread safe.
 * This function doesn't need to be called explicitly.
 */
void random_init();

/**
 * @brief Force another call to srand using current time.
 */
void random_force_reseed();

/**
 * @brief Return an unsigned integer within given range.
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
 * @brief Generate a random polynomial.
 *
 * May also call random_init().
 * poly must be initialized beforehand.
 *
 * @param poly pointer to a polynomial to store the result in
 * @param size maximum size of the polynomial, size must be less or equal to the polynomial->capacity
 * @param weight number of nonzero coefficients, weight <= size
 */
void random_gf4_poly(gf4_poly_t * poly, size_t size);

/**
 * @brief Generate a polynomial of given weight.
 *
 * May also call random_init().
 * poly must be initialized beforehand.
 *
 * @param poly pointer to a polynomial to store the result in
 * @param size maximum size of the polynomial, size <= polynomial->capacity
 * @param weight number of nonzero coefficients, weight <= size
 */
void random_weighted_gf4_poly(gf4_poly_t * poly, size_t size, size_t weight);


/**
 * @brief Generate a polynomial of given weight such that at least weight/2 ones are placed exactly distance apart.
 *
 * poly must be initialized beforehand.
 *
 * @param poly pointer to a polynomial to store the result in
 * @param size maximum size of the polynomial, size <= polynomial->capacity
 * @param weight number of nonzero coefficients, weight <= size
 * @param distance distance between pairs
 */
void random_weighted_gf4_poly_pairs_of_ones(gf4_poly_t * poly, size_t size, size_t weight, size_t distance);

/**
 * @brief Generate a polynomial of given weight such that at least weight/2 pairs of 1 and alpha are placed exactly distance apart.
 *
 * poly must be initialized beforehand.
 *
 * @param poly pointer to a polynomial to store the result in
 * @param size maximum size of the polynomial, size <= polynomial->capacity
 * @param weight number of nonzero coefficients, weight <= size
 * @param distance distance between pairs
 */
void random_weighted_gf4_poly_pairs_of_one_alpha(gf4_poly_t * poly, size_t size, size_t weight, size_t distance);

/**
 * @brief Generate a polynomial of given weight such that at least weight/2 pairs of alpha and 1 are placed exactly distance apart.
 *
 * poly must be initialized beforehand.
 *
 * @param poly pointer to a polynomial to store the result in
 * @param size maximum size of the polynomial, size <= polynomial->capacity
 * @param weight number of nonzero coefficients, weight <= size
 * @param distance distance between pairs
 */
void random_weighted_gf4_poly_pairs_of_alpha_one(gf4_poly_t * poly, size_t size, size_t weight, size_t distance);

#endif //MDPC_GF4_RANDOM_H
