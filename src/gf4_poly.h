#ifndef GF4_GF4_POLY_H
#define GF4_GF4_POLY_H

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "gf4.h"

/**
 * Structure that represents a polynomial over GF4.
 *
 * A polynomial is represented as an array of coefficients.
 * coefficients[0] is the coefficient of x^0.
 * coefficient[degree] is the coefficient of x^degree.
 * It holds that capacity >= degree + 1.
 * In Debug builds, operations over polynomials may resize coefficients array.
 * In Release builds, no resizing can occur! You must initialize polynomials with sufficient capacity in advance!
 */
typedef struct {
	gf4_t * coefficients; ///< an array of coefficients
	size_t degree; ///< degree of the polynomial
	size_t capacity; ///< number of allocated elements in the coefficients array
	                 ///< capacity >= degree+1
} gf4_poly_t;


// initialization
/**
 * Initialize a zero polynomial with the given capacity.
 *
 * Initialized polynomial must be cleaned up using gf4_deinit function if no longer needed!
 *
 * @param capacity number of coefficients to allocate
 * @return initialized polynomial
 */
gf4_poly_t gf4_poly_init_zero(size_t capacity);

/**
 * Zero out an initialized polynomial.
 *
 * Set all coefficients to 0 and set degree to 0. Capacity remains unchanged.
 *
 * @param poly an initialized polynomial
 */
void gf4_poly_zero_out(gf4_poly_t * poly);

/**
 * Destroy a polynomial.
 *
 * Free coefficients array of an initialized polynomial.
 *
 * @param poly an initialized polynomial
 */
void gf4_poly_deinit(gf4_poly_t * poly);

#ifdef SAFE // can resize

/**
 * Change the capacity of an initialized polynomial.
 *
 * It is assumed that new_capacity > poly->capacity. Interally, realloc is called.
 *
 * @param poly an initialized polynomial
 * @param new_capacity required capacity after resizing
 */
void gf4_poly_resize(gf4_poly_t * poly, size_t new_capacity);
#endif // SAFE


// coefficient access
/**
 * Get the coefficient of x^deg.
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
 * Set the coefficient of x^deg to the value val.
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
 * out = a + b
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
 * a = a + b
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
 * poly = poly + val*x^deg
 *
 * poly must be initialized beforehand.
 * The following must hold: deg <= poly->degree.
 *
 * @param poly an initialized polynomial
 * @param deg degree of x
 * @param val coefficient of x^deg
 */
void gf4_poly_add_ax_to_deg_inplace(gf4_poly_t * poly, size_t deg, gf4_t val);


// multiplication
/**
 * out = a * b
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
void gf4_poly_div_x_to_deg(gf4_poly_t * out, gf4_poly_t * poly, size_t deg);
void gf4_poly_div_x_to_deg_inplace(gf4_poly_t * poly, size_t deg);
void gf4_poly_div_rem(gf4_poly_t * div, gf4_poly_t * rem, gf4_poly_t * a, gf4_poly_t * b);

// inverse
bool gf4_poly_invert_slow(gf4_poly_t * maybe_inverse, gf4_poly_t * poly, gf4_poly_t * modulus);

// properties
bool gf4_poly_is_zero(gf4_poly_t * poly);
bool gf4_poly_equal(gf4_poly_t * poly1, gf4_poly_t * poly2);


// helpers
void gf4_poly_adjust_degree(gf4_poly_t * poly, size_t max_degree);
void gf4_poly_pretty_print(gf4_poly_t * poly, FILE * stream, const char * end);
void gf4_poly_coeff_print(gf4_poly_t * poly, size_t max, FILE * stream, const char * end);
gf4_poly_t gf4_poly_clone(gf4_poly_t * poly);
void gf4_poly_copy(gf4_poly_t * out, gf4_poly_t * in);
#endif // GF4_GF4_POLY_H
