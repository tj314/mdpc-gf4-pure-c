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

#include "gf4_poly.h"

// internal
#define GF_4_POLY_GET_DEGREE_INTERNAL(poly) ((poly)->length - 1)

// initialization
gf4_poly_t gf4_poly_init_zero(size_t capacity) {
    assert(1 <= capacity);
    return gf4_vector_init_with_length(capacity, 1); // length=1 => degree=0
}

void gf4_poly_zero_out(gf4_poly_t * poly) {
    assert(NULL != poly);
    memset(poly->array, 0, poly->capacity);
    poly->length = 1;
}

void gf4_poly_deinit(gf4_poly_t * poly) {
    gf4_vector_deinit(poly);
}


// coefficient access
gf4_t gf4_poly_get_coefficient(gf4_poly_t * poly, size_t deg) {
    assert(NULL != poly);
    if (deg > GF_4_POLY_GET_DEGREE_INTERNAL(poly)) {
        return 0;
    } else {
        return poly->array[deg];
    }
}

void gf4_poly_set_coefficient(gf4_poly_t * poly, size_t deg, gf4_t val) {
    assert(NULL != poly);
    if (deg == GF_4_POLY_GET_DEGREE_INTERNAL(poly) && 0 == val) {
        poly->array[deg] = 0;
        gf4_poly_adjust_degree(poly, deg);
    } else if (deg > GF_4_POLY_GET_DEGREE_INTERNAL(poly)) {
        if (0 != val) {
#ifdef CANRESIZE // can resize
            if (deg >= poly->capacity) {
                gf4_vector_resize(poly, deg+1);
			}
#endif // CANRESIZE
            poly->length = deg + 1;
            poly->array[deg] = val;
        }
    } else {
        poly->array[deg] = val;
    }
}


// addition
void gf4_poly_add(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != out);
	assert(NULL != a);
	assert(NULL != b);
    if (GF_4_POLY_GET_DEGREE_INTERNAL(a) == GF_4_POLY_GET_DEGREE_INTERNAL(b)) {
        out->length = 1;
        for (size_t i = 0; i < a->length; ++i) {
            out->array[i] = a->array[i] ^ b->array[i];
            if (0 != out->array[i]) {
                out->length = i+1;
            }
        }
    } else {
        gf4_poly_t * longer;
        gf4_poly_t * shorter;
        if (a->length > b->length) {
            longer = a;
            shorter = b;
        } else {
            longer = b;
            shorter = a;
        }
        out->length = longer->length;
        for (size_t i = 0; i < shorter->length; ++i) {
            out->array[i] = longer->array[i] ^ shorter->array[i];
        }
        for (size_t i = shorter->length; i < longer->length; ++i) {
            out->array[i] = longer->array[i];
        }
    }
}

void gf4_poly_add_inplace(gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != a);
    assert(NULL != b);
    assert(a->capacity >= b->capacity);
    for (size_t i = 0; i < b->length; ++i) {
        a->array[i] ^= b->array[i];
    }
    if (b->length > a->length) {
        a->length = b->length;
    } else if (b->length == a->length) {
        gf4_poly_adjust_degree(a, GF_4_POLY_GET_DEGREE_INTERNAL(a));
    }
}

void gf4_poly_add_ax_to_deg_inplace(gf4_poly_t * poly, size_t deg, gf4_t a) {
    assert(NULL != poly);
#ifdef CANRESIZE
    if (deg >= poly->capacity) {
        gf4_vector_resize(poly, deg + 1);
    }
#endif // CANRESIZE
    poly->array[deg] = poly->array[deg] ^ a;
    if (deg > GF_4_POLY_GET_DEGREE_INTERNAL(poly) && 0 != poly->array[deg]) {
        poly->length = deg + 1;
    } else if (deg == GF_4_POLY_GET_DEGREE_INTERNAL(poly) && 0 == poly->array[deg]) {
        gf4_poly_adjust_degree(poly, deg);
    }
}

// multiplication
void gf4_poly_mul(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != out);
	assert(NULL != a);
	assert(NULL != b);
    for (size_t i = 0; i < a->length; ++i) {
        for (size_t j = 0; j < b->length; ++j) {
            out->array[i+j] ^= gf4_mul(a->array[i], b->array[j]);
        }
    }
    gf4_poly_adjust_degree(out, GF_4_POLY_GET_DEGREE_INTERNAL(a) + GF_4_POLY_GET_DEGREE_INTERNAL(b));
}


// division
void gf4_poly_div_x_to_deg(gf4_poly_t * out, gf4_poly_t * poly, size_t deg) {
    assert(NULL != out);
	assert(NULL != poly);
	assert(GF_4_POLY_GET_DEGREE_INTERNAL(poly) >= deg);
    size_t diff = GF_4_POLY_GET_DEGREE_INTERNAL(poly) - deg;
    for (size_t i = 0; i <= diff; ++i) {
        out->array[i] = poly->array[i + deg];
    }
    for (size_t i = diff + 1; i < poly->length; ++i) {
        out->array[i] = 0;
    }
    out->length = diff + 1;
}

void gf4_poly_div_x_to_deg_inplace(gf4_poly_t * poly, size_t deg) {
	assert(NULL != poly);
	assert(GF_4_POLY_GET_DEGREE_INTERNAL(poly) >= deg);
    size_t diff = GF_4_POLY_GET_DEGREE_INTERNAL(poly) - deg;
    for (size_t i = 0; i <= diff; ++i) {
        poly->array[i] = poly->array[i + deg];
    }
    for (size_t i = diff + 1; i < poly->length; ++i) {
        poly->array[i] = 0;
    }
    poly->length = diff + 1;
}

void gf4_poly_div_rem(gf4_poly_t * div, gf4_poly_t * rem, gf4_poly_t * a, gf4_poly_t * b) {
    assert(!gf4_poly_is_zero(b));

    if (a->length < b->length) {
        gf4_poly_zero_out(div);
        gf4_poly_copy(rem, a);
        return;
    } else {
        gf4_poly_zero_out(div);
        gf4_poly_copy(rem, a);

        gf4_t b_lead = b->array[GF_4_POLY_GET_DEGREE_INTERNAL(b)];

        while (!gf4_poly_is_zero(rem) && GF_4_POLY_GET_DEGREE_INTERNAL(rem) >= GF_4_POLY_GET_DEGREE_INTERNAL(b)) {
            size_t deg = GF_4_POLY_GET_DEGREE_INTERNAL(rem) - GF_4_POLY_GET_DEGREE_INTERNAL(b);
            gf4_t tmp = gf4_div(rem->array[GF_4_POLY_GET_DEGREE_INTERNAL(rem)], b_lead);
            gf4_poly_add_ax_to_deg_inplace(div, deg, tmp);

            for (size_t i = 0; i < b->length; ++i) {
                gf4_t delta = gf4_mul(b->array[i], tmp);
                rem->array[i+deg] ^= delta;
            }

            if (0 == rem->array[GF_4_POLY_GET_DEGREE_INTERNAL(rem)]) {
                gf4_poly_adjust_degree(rem, GF_4_POLY_GET_DEGREE_INTERNAL(rem));
            }
        }
        return;
    }
}


// invert
bool gf4_poly_invert_slow(gf4_poly_t * maybe_inverse, gf4_poly_t * poly, gf4_poly_t * modulus) {
    assert(NULL != maybe_inverse);
    assert(NULL != poly);
    assert(NULL != modulus);
    assert(!gf4_poly_is_zero(modulus));

    if (gf4_poly_is_zero(poly)) {
        return false;
    }

    gf4_poly_t a = gf4_poly_clone(poly);
    gf4_poly_t b = gf4_poly_clone(modulus);
#ifdef CANRESIZE
    gf4_poly_t s = gf4_poly_init_zero(1);
    gf4_poly_t t = gf4_poly_init_zero(1);
    gf4_poly_t div = gf4_poly_init_zero(1);
    gf4_poly_t rem = gf4_poly_init_zero(1);
#else
    gf4_poly_t s = gf4_poly_init_zero(modulus->capacity);
    gf4_poly_t t = gf4_poly_init_zero(modulus->capacity);
    gf4_poly_t div = gf4_poly_init_zero(modulus->capacity);
    gf4_poly_t rem = gf4_poly_init_zero(modulus->capacity);
#endif
    gf4_poly_t tmp;
    s.array[0] = 1;
    while (!gf4_poly_is_zero(&b)) {
        gf4_poly_div_rem(&div, &rem, &a, &b);

        // (a, b) = (b, rem)
        tmp = a;
        a = b;
        b = rem;

        // s = t
        rem = s;
        s = t;

        // t = s + div * t
        gf4_poly_zero_out(&tmp);
        gf4_poly_mul(&tmp, &div, &t);
        gf4_poly_zero_out(&div);
        gf4_poly_add(&div, &rem, &tmp);
        t = div;
        div = tmp;
        gf4_poly_zero_out(&div);
        gf4_poly_zero_out(&rem);
    }

    bool ret_value = false;
    if(0 == GF_4_POLY_GET_DEGREE_INTERNAL(&a)) {
        gf4_t a_abs = a.array[0];
        for (size_t i = 0; i < s.length; ++i) {
            s.array[i] = gf4_div(s.array[i], a_abs);
        }
        gf4_poly_copy(maybe_inverse, &s);
        ret_value = true;
    }
    gf4_poly_deinit(&s);
    gf4_poly_deinit(&a);
    gf4_poly_deinit(&b);
    gf4_poly_deinit(&t);
    gf4_poly_deinit(&div);
    gf4_poly_deinit(&rem);
    return ret_value;
}

// properties
bool gf4_poly_is_zero(gf4_poly_t * poly) {
    assert(NULL != poly);
    return 0 == GF_4_POLY_GET_DEGREE_INTERNAL(poly) && 0 == poly->array[0];
}

bool gf4_poly_equal(gf4_poly_t * poly1, gf4_poly_t * poly2) {
    assert(NULL != poly1);
    assert(NULL != poly2);
    if (poly1->length != poly2->length) {
        return false;
    } else {
        for (size_t i = 0; i < poly1->length; ++i) {
            if (poly1->array[i] != poly2->array[i]) {
                return false;
            }
        }
        return true;
    }
}



// cyclic shifts
void gf4_poly_cyclic_shift_right_inplace(gf4_poly_t * poly, size_t size) {
    assert(NULL != poly);
    assert(size <= poly->capacity);
    gf4_t last = poly->array[size - 1];
    for (size_t i = size - 1; i > 0; --i) {
        poly->array[i] = poly->array[i-1];
    }
    poly->array[0] = last;
}



// helpers
size_t gf4_poly_get_degree(gf4_poly_t * poly) {
    assert(NULL != poly);
    return GF_4_POLY_GET_DEGREE_INTERNAL(poly);
}

void gf4_poly_adjust_degree(gf4_poly_t * poly, size_t max_degree) {
    assert(NULL != poly);
    assert(max_degree < poly->capacity);
    poly->length = 1;
    for (size_t i = max_degree; i != 0; --i) {
        if (0 != poly->array[i]) {
            poly->length = i + 1;
            break;
        }
    }
}

void gf4_poly_pretty_print(gf4_poly_t * poly, FILE * stream, const char * end) {
    assert(NULL != poly);
    assert(NULL != end);
    if (0 == GF_4_POLY_GET_DEGREE_INTERNAL(poly) && 0 == poly->array[0]) {
        fprintf(stream, "0\n");
    } else {
        char * sep = "";
        if (1 <= poly->array[0]) {
            fprintf(stream, "%s", gf4_to_str(poly->array[0]));
            sep = " + ";
        }
        if (1 <= GF_4_POLY_GET_DEGREE_INTERNAL(poly) && 1 == poly->array[1]) {
            fprintf(stream, "%sx", sep);
            sep = " + ";
        } else if (1 <= GF_4_POLY_GET_DEGREE_INTERNAL(poly) && 1 < poly->array[1]) {
            fprintf(stream, "%s%s*x", sep, gf4_to_str(poly->array[1]));
            sep = " + ";
        }
        for (size_t i = 2; i < poly->length; ++i) {
            if (1 == poly->array[i]) {
                fprintf(stream, "%sx^%zu", sep, i);
                sep = " + ";
            } else if (1 < poly->array[i]) {
                fprintf(stream, "%s%s*x^%zu", sep, gf4_to_str(poly->array[i]), i);
                sep = " + ";
            }
        }
        fprintf(stream, "%s", end);
    }
}

gf4_poly_t gf4_poly_clone(gf4_poly_t * poly) {
    assert(NULL != poly);
    gf4_poly_t clone;
    clone.length = poly->length;
    clone.capacity = poly->capacity;
    clone.array = calloc(clone.capacity, sizeof(gf4_t));
    if (NULL == clone.array) {
        fprintf(stderr, "gf4_poly_clone: Memory allocation failed!\n");
        exit(-1);
    }
    memcpy(clone.array, poly->array, sizeof(gf4_t) * clone.length);
    return clone;
}

void gf4_poly_copy(gf4_poly_t * out, gf4_poly_t * in) {
    assert(NULL != out);
    assert(NULL != in);
#ifdef CANRESIZE
    if (out->capacity < in->capacity) {
        gf4_vector_resize(out, in->capacity);
    }
#endif
    out->length = in->length;
    memcpy(out->array, in->array, sizeof(gf4_t)*out->capacity);
}