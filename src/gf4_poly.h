#ifndef GF4_GF4_POLY_H
#define GF4_GF4_POLY_H

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "gf4.h"

typedef struct {
	gf4_t * coefficients;
	size_t degree;
	size_t capacity;
} gf4_poly_t;


// initialization
gf4_poly_t gf4_poly_init_zero(size_t capacity);
void gf4_poly_init_zero_byref(gf4_poly_t * poly, size_t default_capacity);
void gf4_poly_zero_out(gf4_poly_t * poly);
void gf4_poly_deinit(gf4_poly_t * poly);

#ifdef SAFE // can resize
void gf4_poly_resize(gf4_poly_t * poly, size_t new_capacity);
#endif // SAFE


// coefficient access
gf4_t gf4_poly_get_coefficient(gf4_poly_t * poly, size_t deg);
void gf4_poly_set_coefficient(gf4_poly_t * poly, size_t deg, gf4_t val);

// addition
void gf4_poly_add(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b);
void gf4_poly_add_inplace(gf4_poly_t * a, gf4_poly_t * b);
void gf4_poly_add_ax_to_deg_inplace(gf4_poly_t * poly, size_t deg, gf4_t val);


// multiplication
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
void gf4_poly_pretty_print(gf4_poly_t * poly, FILE * stream, const char * end);
void gf4_poly_coeff_print(gf4_poly_t * poly, size_t max, FILE * stream, const char * end);
gf4_poly_t gf4_poly_clone(gf4_poly_t * poly);
void gf4_poly_copy(gf4_poly_t * out, gf4_poly_t * in);
#endif // GF4_GF4_POLY_H
