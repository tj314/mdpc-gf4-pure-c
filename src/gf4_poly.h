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
gf4_poly_t gf4_poly_init_zero(size_t default_capacity);

void gf4_poly_init_zero_byref(gf4_poly_t * poly, size_t default_capacity);


// coefficient access
gf4_t gf4_poly_get_coefficient(gf4_poly_t * poly, size_t deg);

void gf4_poly_set_coefficient(gf4_poly_t * poly, size_t deg, gf4_t val);

// addition
gf4_poly_t gf4_poly_add(gf4_poly_t * a, gf4_poly_t * b);

void gf4_poly_add_byref(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b);


// multiplication
gf4_poly_t gf4_poly_mul(gf4_poly_t * a, gf4_poly_t * b);

void gf4_poly_mul_byref(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b);

// division and modulo
gf4_poly_t gf4_poly_div_x_to_deg(gf4_poly_t * poly, size_t deg);
void gf4_poly_div_x_to_deg_byref(gf4_poly_t * out, gf4_poly_t * poly, size_t deg);
void gf4_poly_div_x_to_deg_inplace(gf4_poly_t * poly, size_t deg);
void gf4_poly_div_rem(gf4_poly_t * div, gf4_poly_t * rem, gf4_poly_t * a, gf4_poly_t * b);

// properties
bool gf4_poly_is_zero(gf4_poly_t * poly);


// helpers
void gf4_poly_pretty_print(gf4_poly_t * poly);
gf4_poly_t gf4_poly_clone(gf4_poly_t * poly);
#endif // GF4_GF4_POLY_H
