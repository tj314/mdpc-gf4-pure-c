#include "gf4_poly.h"

// initialization
gf4_poly_t gf4_poly_init_zero(size_t capacity) {
    gf4_poly_t poly;
    poly.coefficients = calloc(capacity, sizeof(gf4_t));
    if (NULL == poly.coefficients) {
        fprintf(stderr, "gf4_poly_init_zero: Memory allocation failed!\n");
        exit(-1);
    }
    poly.capacity = capacity;
    poly.degree = 0;
    return poly;
}

#ifdef CANRESIZE
void gf4_poly_resize(gf4_poly_t * poly, size_t new_capacity) {
    assert(NULL != poly);
    gf4_t * tmp = realloc(poly->coefficients, new_capacity*sizeof(gf4_t));
    if (NULL == tmp) {
        free(poly->coefficients);
        fprintf(stderr, "gf4_poly_resize: Memory reallocation failed!\n");
        exit(-1);
    }
    poly->coefficients = tmp;
    poly->capacity = new_capacity;
}
#endif

void gf4_poly_zero_out(gf4_poly_t * poly) {
    assert(NULL != poly);
    memset(poly->coefficients, 0, poly->capacity);
    poly->degree = 0;
}

void gf4_poly_deinit(gf4_poly_t * poly) {
    assert(NULL != poly);
    free(poly->coefficients);
    poly->coefficients = NULL;
}


// coefficient access
gf4_t gf4_poly_get_coefficient(gf4_poly_t * poly, size_t deg) {
    assert(NULL != poly);
    if (deg > poly->degree) {
        return 0;
    } else {
        return poly->coefficients[deg];
    }
}

void gf4_poly_set_coefficient(gf4_poly_t * poly, size_t deg, gf4_t val) {
    assert(NULL != poly);
    if (deg == poly->degree && 0 == val) {
        poly->coefficients[deg] = 0;
        gf4_poly_adjust_degree(poly, deg);
    } else if (deg > poly->degree) {
        if (0 != val) {
#ifdef CANRESIZE // can resize
            if (deg >= poly->capacity) {
                gf4_poly_resize(poly, deg+1);
			}
#endif // CANRESIZE
            poly->degree = deg;
            poly->coefficients[deg] = val;
        }
    } else {
        poly->coefficients[deg] = val;
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
            out->coefficients[i] = a->coefficients[i] ^ b->coefficients[i];
            if (0 != out->coefficients[i]) {
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
            out->coefficients[i] = longer->coefficients[i] ^ shorter->coefficients[i];
        }
        for (size_t i = shorter->degree + 1; i <= longer->degree; ++i) {
            out->coefficients[i] = longer->coefficients[i];
        }
    }
}

void gf4_poly_add_inplace(gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != a);
    assert(NULL != b);
    assert(a->capacity >= b->capacity);
    for (size_t i = 0; i <= b->degree; ++i) {
        a->coefficients[i] ^= b->coefficients[i];
    }
    if (b->degree > a->degree) {
        a->degree = b->degree;
    } else if (b->degree == a->degree) {
        gf4_poly_adjust_degree(a, a->degree);
    }
}

void gf4_poly_add_ax_to_deg_inplace(gf4_poly_t * poly, size_t deg, gf4_t val) {
    assert(NULL != poly);
#ifdef CANRESIZE
    if (deg >= poly->capacity) {
        gf4_poly_resize(poly, deg + 1);
    }
#endif // CANRESIZE
    poly->coefficients[deg] = poly->coefficients[deg] ^ val;
    if (deg > poly->degree && 0 != poly->coefficients[deg]) {
        poly->degree = deg;
    } else if (deg == poly->degree && 0 == poly->coefficients[deg]) {
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
            out->coefficients[i+j] ^= gf4_mul(a->coefficients[i], b->coefficients[j]);
        }
    }
    gf4_poly_adjust_degree(out, a->degree + b->degree);
}


// division
void gf4_poly_div_x_to_deg(gf4_poly_t * out, gf4_poly_t * poly, size_t deg) {
    assert(NULL != out);
	assert(NULL != poly);
	assert(poly->degree >= deg);
    size_t diff = poly->degree - deg;
    for (size_t i = 0; i <= diff; ++i) {
        out->coefficients[i] = poly->coefficients[i + deg];
    }
    for (size_t i = diff + 1; i <= poly->degree; ++i) {
        out->coefficients[i] = 0;
    }
    out->degree = diff;
}

void gf4_poly_div_x_to_deg_inplace(gf4_poly_t * poly, size_t deg) {
	assert(NULL != poly);
	assert(poly->degree >= deg);
    size_t diff = poly->degree - deg;
    for (size_t i = 0; i <= diff; ++i) {
        poly->coefficients[i] = poly->coefficients[i + deg];
    }
    for (size_t i = diff + 1; i <= poly->degree; ++i) {
        poly->coefficients[i] = 0;
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

        gf4_t b_lead = b->coefficients[b->degree];

        while (!gf4_poly_is_zero(rem) && rem->degree >= b->degree) {
            size_t deg = rem->degree - b->degree;
            gf4_t tmp = gf4_div(rem->coefficients[rem->degree], b_lead);
            gf4_poly_add_ax_to_deg_inplace(div, deg, tmp);

            for (size_t i = 0; i <= b->degree; ++i) {
                gf4_t delta = gf4_mul(b->coefficients[i], tmp);
                rem->coefficients[i+deg] ^= delta;
            }

            if (0 == rem->coefficients[rem->degree]) {
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
    gf4_poly_t s = gf4_poly_init_zero(modulus->capacity);
    gf4_poly_t t = gf4_poly_init_zero(modulus->capacity);
    gf4_poly_t div = gf4_poly_init_zero(modulus->capacity);
    gf4_poly_t rem = gf4_poly_init_zero(modulus->capacity);
#endif
    gf4_poly_t tmp;
    s.coefficients[0] = 1;
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
        gf4_t a_abs = a.coefficients[0];
        for (size_t i = 0; i <= s.degree; ++i) {
            s.coefficients[i] = gf4_div(s.coefficients[i], a_abs);
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
    return 0 == poly->degree && 0 == poly->coefficients[0];
}

bool gf4_poly_equal(gf4_poly_t * poly1, gf4_poly_t * poly2) {
    assert(NULL != poly1);
    assert(NULL != poly2);
    if (poly1->degree != poly2->degree) {
        return false;
    } else {
        for (size_t i = 0; i <= poly1->degree; ++i) {
            if (poly1->coefficients[i] != poly2->coefficients[i]) {
                return false;
            }
        }
        return true;
    }
}



// helpers
void gf4_poly_adjust_degree(gf4_poly_t * poly, size_t max_degree) {
    assert(NULL != poly);
    poly->degree = 0;
    for (size_t i = max_degree; i != 0; --i) {
        if (0 != poly->coefficients[i]) {
            poly->degree = i;
            break;
        }
    }
}

void gf4_poly_pretty_print(gf4_poly_t * poly, FILE * stream, const char * end) {
    assert(NULL != poly);
    assert(NULL != end);
    if (0 == poly->degree && 0 == poly->coefficients[0]) {
        fprintf(stream, "0\n");
    } else {
        char * sep = "";
        if (1 <= poly->coefficients[0]) {
            fprintf(stream, "%s", gf4_to_str(poly->coefficients[0]));
            sep = " + ";
        }
        if (1 <= poly->degree && 1 == poly->coefficients[1]) {
            fprintf(stream, "%sx", sep);
            sep = " + ";
        } else if (1 <= poly->degree && 1 < poly->coefficients[1]) {
            fprintf(stream, "%s%s*x", sep, gf4_to_str(poly->coefficients[1]));
            sep = " + ";
        }
        for (size_t i = 2; i <= poly->degree; ++i) {
            if (1 == poly->coefficients[i]) {
                fprintf(stream, "%sx^%zu", sep, i);
                sep = " + ";
            } else if (1 < poly->coefficients[i]) {
                fprintf(stream, "%s%s*x^%zu", sep, gf4_to_str(poly->coefficients[i]), i);
                sep = " + ";
            }
        }
        fprintf(stream, "%s", end);
    }
}

void gf4_poly_coeff_print(gf4_poly_t * poly, size_t max, FILE * stream, const char * end) {
    for (size_t i = 0; i < max; ++i) {
        fprintf(stream, "%u ", poly->coefficients[i]);
    }
    fprintf(stream, "%s", end);
}

gf4_poly_t gf4_poly_clone(gf4_poly_t * poly) {
    assert(NULL != poly);
    gf4_poly_t clone;
    clone.degree = poly->degree;
    clone.capacity = poly->capacity;
    clone.coefficients = calloc(clone.capacity, sizeof(gf4_t));
    if (NULL == clone.coefficients) {
        fprintf(stderr, "gf4_poly_clone: Memory allocation failed!\n");
        exit(-1);
    }
    memcpy(clone.coefficients, poly->coefficients, clone.degree+1);
    return clone;
}

void gf4_poly_copy(gf4_poly_t * out, gf4_poly_t * in) {
    assert(NULL != out);
    assert(NULL != in);
#ifdef CANRESIZE
    if (out->capacity < in->capacity) {
        gf4_poly_resize(out, in->capacity);
    }
#endif
    out->degree = in->degree;
    memcpy(out->coefficients, in->coefficients, out->capacity);
}