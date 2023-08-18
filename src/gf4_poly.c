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
#include "utils.h"

// initialization
gf4_poly_t gf4_poly_init_zero(size_t capacity) {
    assert(1 <= capacity);
    gf4_poly_t poly;
    poly.coefficients = gf4_array_init(capacity, true);
    poly.degree = 0;
    return poly;
}

void gf4_poly_zero_out(gf4_poly_t * poly) {
    assert(NULL != poly);
    memset(poly->coefficients.array, 0, sizeof(gf4_t)*poly->coefficients.capacity);
    poly->degree = 0;
}

void gf4_poly_deinit(gf4_poly_t * poly) {
    gf4_array_deinit(&poly->coefficients);
    poly->degree = 0;
}


// coefficient access
gf4_t gf4_poly_get_coefficient(gf4_poly_t * poly, size_t deg) {
    assert(NULL != poly);
    if (deg > poly->degree) {
        return 0;
    } else {
        return poly->coefficients.array[deg];
    }
}

void gf4_poly_set_coefficient(gf4_poly_t * poly, size_t deg, gf4_t val) {
    assert(NULL != poly);
    if (deg == poly->degree && 0 == val) {
        poly->coefficients.array[deg] = 0;
        gf4_poly_adjust_degree(poly, deg);
    } else if (deg > poly->degree) {
        if (0 != val) {
#ifdef CANRESIZE // can resize
            if (deg >= poly->coefficients.capacity) {
                gf4_array_resize(&poly->coefficients, deg+1);
			}
#endif // CANRESIZE
            poly->degree = deg;
            poly->coefficients.array[deg] = val;
        }
    } else {
        poly->coefficients.array[deg] = val;
    }
}


// addition
void gf4_poly_add(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != out);
	assert(NULL != a);
	assert(NULL != b);
    if (a->degree == b->degree) {
        out->degree = 0;
        for (size_t i = 0; i <= a->degree; ++i) {
            out->coefficients.array[i] = a->coefficients.array[i] ^ b->coefficients.array[i];
            if (0 != out->coefficients.array[i]) {
                out->degree = i;
            }
        }
    } else {
        gf4_poly_t * longer;
        gf4_poly_t * shorter;
        if (a->degree > b->degree) {
            longer = a;
            shorter = b;
        } else {
            longer = b;
            shorter = a;
        }
        out->degree = longer->degree;
        for (size_t i = 0; i <= shorter->degree; ++i) {
            out->coefficients.array[i] = longer->coefficients.array[i] ^ shorter->coefficients.array[i];
        }
        for (size_t i = shorter->degree + 1; i <= longer->degree; ++i) {
            out->coefficients.array[i] = longer->coefficients.array[i];
        }
    }
}

void gf4_poly_add_inplace(gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != a);
    assert(NULL != b);
    assert(a->coefficients.capacity >= b->coefficients.capacity);
    for (size_t i = 0; i <= b->degree; ++i) {
        a->coefficients.array[i] ^= b->coefficients.array[i];
    }
    if (b->degree > a->degree) {
        a->degree = b->degree;
    } else if (b->degree == a->degree) {
        gf4_poly_adjust_degree(a, a->degree);
    }
}

void gf4_poly_add_ax_to_deg_inplace(gf4_poly_t * poly, size_t deg, gf4_t a) {
    assert(NULL != poly);
#ifdef CANRESIZE
    if (deg >= poly->coefficients.capacity) {
        gf4_array_resize(&poly.coefficients, deg + 1);
    }
#endif // CANRESIZE
    poly->coefficients.array[deg] = poly->coefficients.array[deg] ^ a;
    if (deg > poly->degree && 0 != poly->coefficients.array[deg]) {
        poly->degree = deg;
    } else if (deg == poly->degree && 0 == poly->coefficients.array[deg]) {
        gf4_poly_adjust_degree(poly, deg);
    }
}

// multiplication
void gf4_poly_mul(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != out);
	assert(NULL != a);
	assert(NULL != b);
    for (size_t i = 0; i <= a->degree; ++i) {
        for (size_t j = 0; j <= b->degree; ++j) {
            out->coefficients.array[i+j] ^= gf4_mul(a->coefficients.array[i], b->coefficients.array[j]);
        }
    }
    gf4_poly_adjust_degree(out, a->degree + b->degree);
}


// division
void gf4_poly_div_x_to_deg(gf4_poly_t * out, gf4_poly_t * poly, size_t deg) {
    assert(NULL != out);
	assert(NULL != poly);
	assert(out->coefficients.capacity > UTILS_SUBTRACT_OR_ZERO(poly->degree, deg));
    size_t diff = UTILS_SUBTRACT_OR_ZERO(poly->degree, deg);
    for (size_t i = 0; i <= diff; ++i) {
        out->coefficients.array[i] = poly->coefficients.array[i + deg];
    }
    for (size_t i = diff + 1; i <= poly->degree; ++i) {
        out->coefficients.array[i] = 0;
    }
    out->degree = diff;
}

void gf4_poly_div_x_to_deg_inplace(gf4_poly_t * poly, size_t deg) {
	assert(NULL != poly);
    assert(poly->coefficients.capacity > UTILS_SUBTRACT_OR_ZERO(poly->degree, deg));
    size_t diff = UTILS_SUBTRACT_OR_ZERO(poly->degree, deg);
    for (size_t i = 0; i <= diff; ++i) {
        poly->coefficients.array[i] = poly->coefficients.array[i + deg];
    }
    for (size_t i = diff + 1; i <= poly->degree; ++i) {
        poly->coefficients.array[i] = 0;
    }
    poly->degree = diff;
}

void gf4_poly_div_rem(gf4_poly_t * div, gf4_poly_t * rem, gf4_poly_t * a, gf4_poly_t * b) {
    assert(!gf4_poly_is_zero(b));

    if (a->degree < b->degree) {
        gf4_poly_zero_out(div);
        gf4_poly_copy(rem, a);
        return;
    } else {
        gf4_poly_zero_out(div);
        gf4_poly_copy(rem, a);

        gf4_t b_lead = b->coefficients.array[b->degree];

        while (!gf4_poly_is_zero(rem) && rem->degree >= b->degree) {
            size_t deg = rem->degree - b->degree;
            gf4_t tmp = gf4_div(rem->coefficients.array[rem->degree], b_lead);
            gf4_poly_add_ax_to_deg_inplace(div, deg, tmp);

            for (size_t i = 0; i <= b->degree; ++i) {
                gf4_t delta = gf4_mul(b->coefficients.array[i], tmp);
                rem->coefficients.array[i+deg] ^= delta;
            }

            if (0 == rem->coefficients.array[rem->degree]) {
                gf4_poly_adjust_degree(rem, rem->degree);
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
    gf4_poly_t s = gf4_poly_init_zero(modulus->coefficients.capacity);
    gf4_poly_t t = gf4_poly_init_zero(modulus->coefficients.capacity);
    gf4_poly_t div = gf4_poly_init_zero(modulus->coefficients.capacity);
    gf4_poly_t rem = gf4_poly_init_zero(modulus->coefficients.capacity);
#endif
    gf4_poly_t tmp;
    s.coefficients.array[0] = 1;
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
    if(0 == a.degree) {
        gf4_t a_abs = a.coefficients.array[0];
        for (size_t i = 0; i <= s.degree; ++i) {
            s.coefficients.array[i] = gf4_div(s.coefficients.array[i], a_abs);
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
    return 0 == poly->degree && 0 == poly->coefficients.array[0];
}

bool gf4_poly_equal(gf4_poly_t * poly1, gf4_poly_t * poly2) {
    assert(NULL != poly1);
    assert(NULL != poly2);
    if (poly1->degree != poly2->degree) {
        return false;
    } else {
        for (size_t i = 0; i <= poly1->degree; ++i) {
            if (poly1->coefficients.array[i] != poly2->coefficients.array[i]) {
                return false;
            }
        }
        return true;
    }
}



// cyclic shifts
void gf4_poly_cyclic_shift_right_inplace(gf4_poly_t * poly, size_t size) {
    assert(NULL != poly);
    assert(size <= poly->coefficients.capacity);
    gf4_t last = poly->coefficients.array[size - 1];
    for (size_t i = size - 1; i > 0; --i) {
        poly->coefficients.array[i] = poly->coefficients.array[i-1];
    }
    poly->coefficients.array[0] = last;
}



// helpers
size_t gf4_poly_get_degree(gf4_poly_t * poly) {
    assert(NULL != poly);
    return poly->degree;
}

void gf4_poly_adjust_degree(gf4_poly_t * poly, size_t max_degree) {
    assert(NULL != poly);
    assert(max_degree < poly->coefficients.capacity);
    poly->degree = 0;
    for (size_t i = max_degree; i != 0; --i) {
        if (0 != poly->coefficients.array[i]) {
            poly->degree = i;
            break;
        }
    }
}

void gf4_poly_pretty_print(gf4_poly_t * poly, FILE * stream, const char * end) {
    assert(NULL != poly);
    assert(NULL != end);
    if (0 == poly->degree && 0 == poly->coefficients.array[0]) {
        fprintf(stream, "0\n");
    } else {
        char * sep = "";
        if (1 <= poly->coefficients.array[0]) {
            fprintf(stream, "%s", gf4_to_str(poly->coefficients.array[0]));
            sep = " + ";
        }
        if (1 <= poly->degree && 1 == poly->coefficients.array[1]) {
            fprintf(stream, "%sx", sep);
            sep = " + ";
        } else if (1 <= poly->degree && 1 < poly->coefficients.array[1]) {
            fprintf(stream, "%s%s*x", sep, gf4_to_str(poly->coefficients.array[1]));
            sep = " + ";
        }
        for (size_t i = 2; i <= poly->degree; ++i) {
            if (1 == poly->coefficients.array[i]) {
                fprintf(stream, "%sx^%zu", sep, i);
                sep = " + ";
            } else if (1 < poly->coefficients.array[i]) {
                fprintf(stream, "%s%s*x^%zu", sep, gf4_to_str(poly->coefficients.array[i]), i);
                sep = " + ";
            }
        }
        fprintf(stream, "%s", end);
    }
}

gf4_poly_t gf4_poly_clone(gf4_poly_t * poly) {
    assert(NULL != poly);
    gf4_poly_t clone;
    clone.degree = poly->degree;
    clone.coefficients = gf4_array_clone(&poly->coefficients);
    return clone;
}

void gf4_poly_copy(gf4_poly_t * out, gf4_poly_t * in) {
    assert(NULL != out);
    assert(NULL != in);
#ifdef CANRESIZE
    if (out->coefficients.capacity < in->coefficients.capacity) {
        gf4_array_resize(&out->coefficients, in->coefficients.capacity);
    }
#endif
    out->degree = in->degree;
    memcpy(out->coefficients.array, in->coefficients.array, sizeof(gf4_t)*out->coefficients.capacity);
}