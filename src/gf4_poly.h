/**
 *  @file   gf4_poly.h
 *  @brief  Polynomials over GF(4).
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

#ifndef GF4_GF4_POLY_H
#define GF4_GF4_POLY_H

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "gf4.h"
#include "gf4_vector.h"

/**
 * @brief Structure that represents a polynomial over GF4.
 *
 * A polynomial is represented as an array of coefficients.
 * array[0] is the coefficient of x^0.
 * The degree of polynomial is length-1.
 * array[length - 1] is the coefficient of x^degree.
 * In Debug builds, operations over polynomials may resize coefficients array.
 * In Release builds, no resizing can occur! You must initialize polynomials with sufficient capacity in advance!
 */
 typedef gf4_vector_t gf4_poly_t;
// initialization
/**
 * @brief Initialize a zero polynomial with the given capacity.
 *
 * Initialized polynomial must be cleaned up using gf4_poly_deinit function if no longer needed!
 *
 * @see gf4_poly_deinit
 *
 * @param capacity number of coefficients to allocate
 * @return initialized polynomial
 */
gf4_poly_t gf4_poly_init_zero(size_t capacity);

/**
 * @brief Zero out an initialized polynomial.
 *
 * Set all coefficients to 0 and set degree to 0. Capacity remains unchanged.
 *
 * @see gf4_poly_deinit
 * @param poly an initialized polynomial
 */
void gf4_poly_zero_out(gf4_poly_t * poly);

/**
 * @brief Destroy a polynomial.
 *
 * Free coefficients array of an initialized polynomial.
 *
 * @param poly an initialized polynomial
 */
void gf4_poly_deinit(gf4_poly_t * poly);

// coefficient access
/**
 * @brief Get the coefficient of x^deg.
 *
 * Equivalent to poly->coefficients[deg] but with added bonus of bounds checking.
 * If deg > poly->degree, 0 is returned.
 * If bounds checking is not desired, it is recommended to use direct access to the coefficients array instead.
 *
 * @param poly an initialized polynomial
 * @param deg degree of x
 * @return coefficient of x^deg or 0 if deg > poly->degree
 */
gf4_t gf4_poly_get_coefficient(gf4_poly_t * poly, size_t deg);

/**
 * @brief Set the coefficient of x^deg to the value val.
 *
 * If resizing is enabled and deg >= capacity, the polynomial is resized.
 * If resizing is not enabled and deg >= capacity, the behavior is undefined.
 * The degree of poly is updated to the correct value after the change.
 *
 * @param poly an initialized polynomial
 * @param deg degree of x
 * @param val value to be set
 */
void gf4_poly_set_coefficient(gf4_poly_t * poly, size_t deg, gf4_t val);

// addition

/**
 * @brief out = a + b
 *
 * out, a and b must be initialized beforehand.
 * If resizing is enabled, out may be resized to fit the result of addition.
 * If it is not enabled, out must have capacity > max(a->degree, b->degree).
 *
 * @param out pointer to a polynomial to store the result in
 * @param a pointer to a polynomial
 * @param b pointer to a polynomial
 */
void gf4_poly_add(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b);

/**
 * @brief a = a + b
 *
 * a and b must be initialized beforehand.
 * If resizing is enabled, a may be resized to fit the result of addition.
 * If it is not enabled, a must have capacity > b->degree.
 *
 * @param a pointer to a polynomial
 * @param b pointer to a polynomial
 */
void gf4_poly_add_inplace(gf4_poly_t * a, gf4_poly_t * b);

/**
 * @brief poly = poly + a*x^deg
 *
 * poly must be initialized beforehand.
 * The following must hold: deg <= poly->degree.
 *
 * @param poly an initialized polynomial
 * @param deg degree of x
 * @param a coefficient of x^deg
 */
void gf4_poly_add_ax_to_deg_inplace(gf4_poly_t * poly, size_t deg, gf4_t a);


// multiplication
/**
 * @brief out = a * b
 *
 * out, a and b must be initialized beforehand.
 * If resizing is enabled, out may be resized to fit the result of multiplication.
 * If it is not enabled, out must have capacity > a->degree + b->degree.
 *
 * @param out pointer to a polynomial to store the result in
 * @param a pointer to a polynomial
 * @param b pointer to a polynomial
 */
void gf4_poly_mul(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b);

// division and modulo
/**
 * @brief out = poly / x^deg
 *
 * out and poly must be initialized beforehand.
 *
 * If resizing is enabled, out may be resized to fit the result of multiplication.
 * If it is not enabled, out must have capacity > max(poly->degree - deg, 0).
 *
 * @param out pointer to a polynomial to be divided
 * @param poly polynomial to be divided
 * @param deg degree of x^deg
 */
void gf4_poly_div_x_to_deg(gf4_poly_t * out, gf4_poly_t * poly, size_t deg);

/**
 * @brief poly = poly / x^deg
 *
 * poly must be initialized beforehand.
 *
 * @param poly polynomial to be divided
 * @param deg degree of x^deg
 */
void gf4_poly_div_x_to_deg_inplace(gf4_poly_t * poly, size_t deg);

/**
 * @brief (div, mod) = (a / b, a % b)
 *
 * div, rem, a, b must be initialized beforehand.
 * If resizing is enabled, div and rem may be resized to fit the result of multiplication.
 * If it is not enabled, div and rem must have sufficient capacity.
 *
 * @param div pointer to a polynomial to store the result of division in
 * @param rem pointer to a polynomial to store the resulting remainder in
 * @param a pointer to a polynomial to be used as the dividend
 * @param b pointer to a polynomial to be used as the divisor
 */
void gf4_poly_div_rem(gf4_poly_t * div, gf4_poly_t * rem, gf4_poly_t * a, gf4_poly_t * b);

// inverse
/**
 * @brief maybe_inverse = poly^-1 (mod modulus)
 *
 * Implemented using xgcd.
 * maybe_inverse will be zeroed out if the inverse does not exist.
 *
 * maybe_inverse, poly, modulus must be initialized beforehand.
 * If resizing is enabled, maybe_inverse can be resized to fit the inverse.
 * If it is not enabled, maybe_inverse must have sufficient capacity.
 *
 * @param maybe_inverse pointer to the polynomial to store the inverse in
 * @param poly pointer to the polynomial to be inverted
 * @param modulus pointer to the polynomial to by used as a modulus
 * @return true if the inverse was found, false otherwise
 */
bool gf4_poly_invert_slow(gf4_poly_t * maybe_inverse, gf4_poly_t * poly, gf4_poly_t * modulus);

// properties
/**
 * @brief check whether polynomial is zero.
 *
 * poly must be initialized beforehand.
 *
 * @param poly pointer to a polynomial to be checked
 * @return true if poly is a zero polynomial, false otherwise
 */
bool gf4_poly_is_zero(gf4_poly_t * poly);

/**
 * @brief Compare poly1 and poly2.
 *
 * poly1 and poly2 must be initialized beforehand.
 * poly1 and poly2 must have the same capacity.
 *
 * @param poly1 pointer to a polynomial
 * @param poly2 pointer to a polynomial
 * @return
 */
bool gf4_poly_equal(gf4_poly_t * poly1, gf4_poly_t * poly2);

// cyclic shifts
/**
 * @brief Perform a cyclic shift of the coefficients of poly to right by one place.
 *
 * size determines where to overflow coefficients to the beginning of the array,
 * i.e. size is the index of the first element not to be shifted.
 * The following must be true: size <= poly->capacity.
 * e.g.: cyclic shift of 1 + x + x^3 with size=4 and capacity=5 is performed as follows: [1, 1, 0, 1, 0] -> [1, 1, 1, 0, 0] -> 1 + x + x^2
 * e.g.: cyclic shift of 1 + x + x^3 with size=5 and capacity=5 is performed as follows: [1, 1, 0, 1, 0] -> [0, 1, 1, 0, 1] -> x + x^2 + x^4
 *
 * @param poly pointer to a polynomial to be shifted
 * @param size size of the shift, size <= poly->capacity
 */
void gf4_poly_cyclic_shift_right_inplace(gf4_poly_t * poly, size_t size);


// helpers
/**
 * @brief Get the degree of the polynomial.
 *
 * @param poly a pointer to a polynomial
 * @return the degree of this polynomial
 */
size_t gf4_poly_get_degree(gf4_poly_t * poly);

/**
 * @brief Correctly set the degree of a polynomial.
 *
 * poly must be initialized beforehand.
 * If the previous valid degree is known and it is guaranteed that the new degree will be lower,
 * you can set max_degree to it. Otherwise max_degree must be set to poly->capacity - 1.
 * The following must hold: max_degree < poly.capacity.
 *
 * The search for the correct degree must is performed by a loop starting at max_degree and iterating towards 0.
 * The loop breaks if a nonzero coefficient is found.
 *
 * @param poly pointer to a polynomial whose degree is to be adjusted
 * @param max_degree maximum possible value that the new degree could equal
 */
void gf4_poly_adjust_degree(gf4_poly_t * poly, size_t max_degree);

/**
 * @brief Print nicely formatted polynomial.
 *
 * poly must be initialized beforehand.
 * e.g.: [0, 1, 0, 2, 2, 1, 0] --> x + (a+1)*x^3 + (a+1)*x^4 + x^5
 *
 * @param poly pointer to a polynomial to be printed.
 * @param stream stream to be used (e.g. stdout, stderr...)
 * @param end last character to be printed (e. g. line ending or space)
 */
void gf4_poly_pretty_print(gf4_poly_t * poly, FILE * stream, const char * end);

/**
 * @brief Clone a polynomial.
 *
 * poly must be initialized beforehand.
 * A new polynomial is created and returned.
 * This polynomial contains a copy of coefficients of poly. capacity and degree are the same.
 * Cloned polynomial must be cleaned up using gf4_deinit function if no longer needed!
 *
 * @see gf4_poly_deinit
 *
 * @param poly pointer to a polynomial
 * @return copied polynomial
 */
gf4_poly_t gf4_poly_clone(gf4_poly_t * poly);

/**
 * @brief Copy a polynomial.
 *
 * Copy in to out.
 *
 * poly must be initialized beforehand.
 * If resizing is enabled, out can be resized to fit coefficients of in.
 * If it is not enabled, out must have sufficient capacity.
 *
 * @param out pointer to a polynomial
 * @param in pointer to a polynomial
 */
void gf4_poly_copy(gf4_poly_t * out, gf4_poly_t * in);
#endif // GF4_GF4_POLY_H
