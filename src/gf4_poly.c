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

void gf4_poly_init_zero_byref(gf4_poly_t * poly, size_t default_capacity) {
    poly->coefficients = calloc(default_capacity, sizeof(gf4_t));
    if (NULL == poly->coefficients) {
        fprintf(stderr, "gf4_poly_init_byref: Memory allocation failed!\n");
        exit(-1);
    }
    poly->capacity = default_capacity;
    poly->degree = 0;
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
    memset(poly->coefficients, 0, poly->capacity);
}

void gf4_poly_deinit(gf4_poly_t * poly) {
    free(poly->coefficients);
    poly->coefficients = NULL;
}


// coefficient access
gf4_t gf4_poly_get_coefficient(gf4_poly_t * poly, size_t deg) {
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
        for (size_t i = deg; i > 0; --i) {
            if (0 != poly->coefficients[i]) {
                poly->degree = i;
                return;
            }
        }
        poly->degree = 0;
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
gf4_poly_t gf4_poly_add(gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != a);
	assert(NULL != b);
    if (a->degree == b->degree) {
        gf4_poly_t out;
        out.degree = 0;
        out.capacity = a->capacity;
        out.coefficients = calloc(out.capacity, sizeof(gf4_t));
        if (NULL == out.coefficients) {
            fprintf(stderr, "gf4_poly_add: Memory allocation failed!\n");
            exit(-1);
        }
        for (size_t i = 0; i <= a->degree; ++i) {
            out.coefficients[i] = a->coefficients[i] ^ b->coefficients[i];
            if (0 != out.coefficients[i]) {
                out.degree = i;
            }
        }
        return out;
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
        gf4_poly_t out;
        out.degree = longer->degree;
        out.capacity = longer->capacity;
        out.coefficients = calloc(out.capacity, sizeof(gf4_t));
        assert(NULL != out.coefficients);
        for (size_t i = 0; i <= shorter->degree; ++i) {
            out.coefficients[i] = longer->coefficients[i] ^ shorter->coefficients[i];
        }
        for (size_t i = shorter->degree + 1; i <= longer->degree; ++i) {
            out.coefficients[i] = longer->coefficients[i];
        }
        return out;
    }
}

void gf4_poly_add_byref(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
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
        for (size_t i = deg; i > 0; --i) {
            if (0 != poly->coefficients[i]) {
                poly->degree = i;
                return;
            }
        }
        poly->degree = 0;
    }
}

// multiplication
gf4_poly_t gf4_poly_mul(gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != a);
	assert(NULL != b);
    gf4_poly_t out;
    out.capacity = a->degree + b->degree + 1;
    out.capacity = out.capacity > a->capacity ? out.capacity : a->capacity;
    out.capacity = out.capacity > b->capacity ? out.capacity : b->capacity;
    out.coefficients = calloc(out.capacity, sizeof(gf4_t));
    for (size_t i = 0; i <= a->degree; ++i) {
        for (size_t j = 0; j <= b->degree; ++j) {
            out.coefficients[i+j] ^= gf4_mul(a->coefficients[i], b->coefficients[j]);
        }
    }
    out.degree = 0;
    for (size_t i = 0; i < out.capacity; ++i) {
        if (0 != out.coefficients[i]) {
            out.degree = i;
        }
    }
    return out;
}

void gf4_poly_mul_byref(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
    assert(NULL != out);
	assert(NULL != a);
	assert(NULL != b);
    for (size_t i = 0; i <= a->degree; ++i) {
        for (size_t j = 0; j <= b->degree; ++j) {
            out->coefficients[i+j] ^= gf4_mul(a->coefficients[i], b->coefficients[j]);
        }
    }
    out->degree = 0;
    for (size_t i = 0; i < out->capacity; ++i) {
        if (0 != out->coefficients[i]) {
            out->degree = i;
        }
    }
}

// division and modulo
gf4_poly_t gf4_poly_div_x_to_deg(gf4_poly_t * poly, size_t deg) {
    assert(NULL != poly);
	assert(poly->degree >= deg);
    gf4_poly_t out;
    out.capacity = poly->capacity;
    out.coefficients = calloc(out.capacity, sizeof(gf4_t));
    if (NULL == out.coefficients) {
        fprintf(stderr, "gf4_poly_div_x_to_deg: Memory allocation failed!\n");
        exit(-1);
    }

    size_t diff = poly->degree - deg;
    for (size_t i = 0; i <= diff; ++i) {
        out.coefficients[i] = poly->coefficients[i + deg];
    }
    for (size_t i = diff + 1; i <= poly->degree; ++i) {
        out.coefficients[i] = 0;
    }
    out.degree = diff;
    return out;
}

void gf4_poly_div_x_to_deg_byref(gf4_poly_t * out, gf4_poly_t * poly, size_t deg) {
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
#ifdef CANRESIZE
        gf4_poly_init_zero_byref(div, 1);
#else
        gf4_poly_init_zero_byref(div, a->capacity);
#endif // CANRESIZE
        *rem = gf4_poly_clone(a);
        return;
    } else {
        gf4_poly_t q, r;
#ifdef CANRESIZE
        gf4_poly_init_zero_byref(&q, 1);
#else
        gf4_poly_init_zero_byref(&q, a->capacity);
#endif
        r = gf4_poly_clone(a);

        gf4_t b_lead = b->coefficients[b->degree];

        while (!gf4_poly_is_zero(&r) && r.degree >= b->degree) {
            size_t deg = r.degree - b->degree;
            gf4_t tmp = gf4_div(r.coefficients[r.degree], b_lead);
            gf4_poly_add_ax_to_deg_inplace(&q, deg, tmp);

            for (size_t i = 0; i <= b->degree; ++i) {
                gf4_t delta = gf4_mul(b->coefficients[i], tmp);
                r.coefficients[i+deg] ^= delta;
            }

            if (0 == r.coefficients[r.degree]) {
                size_t old_deg = r.degree;
                r.degree = 0;
                for (size_t i = old_deg; i > 0; --i) {
                    if (0 != r.coefficients[i]) {
                        r.degree = i;
                        break;
                    }
                }
            }
        }
        *div = q;
        *rem = r;
        return;
    }
}


// invert
bool gf4_poly_invert_slow_byref(gf4_poly_t * maybe_inverse, gf4_poly_t * poly, gf4_poly_t * modulus) {
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
        gf4_poly_shallow_copy(&tmp, &a);
        gf4_poly_shallow_copy(&a, &b);
        gf4_poly_shallow_copy(&b, &rem);
        // s = t
        gf4_poly_zero_out(&tmp);
        gf4_poly_shallow_copy(&rem, &s);
        gf4_poly_shallow_copy(&s, &t);
        // t = s + div * t
        gf4_poly_mul_byref(&tmp, &div, &t);
        gf4_poly_zero_out(&div);
        gf4_poly_add_byref(&div, &rem, &tmp);
        gf4_poly_shallow_copy(&t, &div);
        gf4_poly_shallow_copy(&div, &tmp);
        gf4_poly_zero_out(&div);
        gf4_poly_zero_out(&rem);
    }

    bool ret_value = false;
    if(a.degree == 0) {
        gf4_t a_abs = a.coefficients[0];
        for (size_t i = 0; i <= s.degree; ++i) {
            s.coefficients[i] = gf4_div(s.coefficients[i], a_abs);
        }
        *maybe_inverse = s;
        ret_value = true;
    }
    if (!ret_value) {
        gf4_poly_deinit(&s);
    }
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
void gf4_poly_pretty_print(gf4_poly_t * poly, const char * end) {
    assert(NULL != poly);
    assert(NULL != end);
    if (0 == poly->degree && 0 == poly->coefficients[0]) {
        printf("0\n");
    } else {
        char * sep = "";
        if (1 == poly->coefficients[0]) {
            printf("1");
            sep = " + ";
        } else if (1 < poly->coefficients[0]) {
            printf("%s", gf4_to_str(poly->coefficients[0]));
            sep = " + ";
        }
        if (1 <= poly->degree && 1 == poly->coefficients[1]) {
            printf("%s1x", sep);
            sep = " + ";
        } else if (1 <= poly->degree && 1 < poly->coefficients[1]) {
            printf("%s%sx", sep, gf4_to_str(poly->coefficients[1]));
            sep = " + ";
        }
        for (size_t i = 2; i <= poly->degree; ++i) {
            if (1 == poly->coefficients[i]) {
                printf("%sx^%zu", sep, i);
                sep = " + ";
            } else if (1 < poly->coefficients[i]) {
                printf("%s%sx^%zu", sep, gf4_to_str(poly->coefficients[i]), i);
                sep = " + ";
            }
        }
        printf("%s", end);
    }
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

void gf4_poly_shallow_copy(gf4_poly_t * out, gf4_poly_t * in) {
    out->degree = in->degree;
    out->capacity = in->capacity;
    out->coefficients = in->coefficients;
}