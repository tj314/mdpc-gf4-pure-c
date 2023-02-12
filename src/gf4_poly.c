#include "gf4_poly.h"

// initialization
gf4_poly_t gf4_poly_init_zero(size_t default_capacity) {
    gf4_poly_t poly;
    poly.coefficients = calloc(default_capacity, sizeof(gf4_t));
    assert(NULL != poly.coefficients);
    poly.capacity = default_capacity;
    poly.degree = 0;
    return poly;
}

void gf4_poly_init_zero_byref(gf4_poly_t * poly, size_t default_capacity) {
    poly->coefficients = calloc(default_capacity, sizeof(gf4_t));
    assert(NULL != poly->coefficients);
    poly->capacity = default_capacity;
    poly->degree = 0;
}

#ifdef SAFE // can resize
void gf4_poly_resize(gf4_poly_t * poly, size_t new_capacity) {
    assert(NULL != poly);
    gf4_t * tmp = realloc(poly->coefficients, new_capacity*sizeof(gf4_t));
    assert(NULL != tmp);
    poly->coefficients = tmp;
    poly->capacity = new_capacity;
}
#endif // SAFE

void gf4_poly_copy_from(gf4_poly_t * out, gf4_poly_t * in) {
#ifdef SAFE
    assert(NULL != out);
    assert(NULL != in);
    if (out->capacity < in->capacity) {
        gf4_poly_resize(out, in->capacity);
    }
#endif
    for (size_t i = 0; i < in->capacity; ++i) {
        out->coefficients[i] = in->coefficients[i];
    }
    for (size_t i = in->capacity; i < out->capacity; ++i) {
        out->coefficients[i] = 0;
    }
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
#ifdef SAFE
    assert(NULL != poly);
#endif // SAFE
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
#ifdef SAFE // can resize
            if (deg >= poly->capacity) {
				gf4_t * tmp = realloc(poly->coefficients, (deg+1)*sizeof(gf4_t));
				assert(NULL != tmp);
				poly->coefficients = tmp;
				poly->capacity = deg+1;
			}
#endif // SAFE
            poly->degree = deg;
            poly->coefficients[deg] = val;
        }
    } else {
        poly->coefficients[deg] = val;
    }
}

// addition
gf4_poly_t gf4_poly_add(gf4_poly_t * a, gf4_poly_t * b) {
#ifdef SAFE
    assert(NULL != a);
	assert(NULL != b);
#endif // SAFE
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
    for (size_t i = shorter->degree; i <= longer->degree; ++i) {
        out.coefficients[i] = longer->coefficients[i];
    }
    return out;
}

void gf4_poly_add_byref(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
#ifdef SAFE
    assert(NULL != out);
	assert(NULL != a);
	assert(NULL != b);
#endif // SAFE
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
    for (size_t i = shorter->degree; i <= longer->degree; ++i) {
        out->coefficients[i] = longer->coefficients[i];
    }
}

void gf4_poly_add_ax_to_deg_inplace(gf4_poly_t * poly, gf4_t a, size_t deg) {
#ifdef SAFE
    assert(NULL != poly);
    gf4_poly_resize(poly, deg + 1);
#endif // SAFE
    poly->coefficients[deg] = poly->coefficients[deg] ^ a;
    if (0 == poly->coefficients[deg]) {
        if (deg > poly->degree) {
            poly->degree = deg;
        } else if (deg == poly->degree) {
            for (size_t i = deg; i > 0; --i) {
                if (0 != poly->coefficients[i]) {
                    poly->degree = i;
                    return;
                }
            }
            poly->degree = 0;
        }
    }
}

// multiplication
gf4_poly_t gf4_poly_mul(gf4_poly_t * a, gf4_poly_t * b) {
#ifdef SAFE
    assert(NULL != a);
	assert(NULL != b);
#endif // SAFE
    gf4_poly_t out;
    out.capacity = a->capacity > b->capacity ? a->capacity : b->capacity;
    out.coefficients = calloc(out.capacity, sizeof(gf4_t));
    for (size_t i = 0; i <= a->degree; ++i) {
        for (size_t j = 0; j <= b->degree; ++j) {
            out.coefficients[i+j] = out.coefficients[i+j] ^ gf4_mul(a->coefficients[i], b->coefficients[j]);
        }
    }
    for (size_t i = 0; i < out.capacity; ++i) {
        if (0 != out.coefficients[i]) {
            out.degree = i;
        }
    }
    return out;
}

void gf4_poly_mul_byref(gf4_poly_t * out, gf4_poly_t * a, gf4_poly_t * b) {
#ifdef SAFE
    assert(NULL != out);
	assert(NULL != a);
	assert(NULL != b);
#endif // SAFE
    for (size_t i = 0; i <= a->degree; ++i) {
        for (size_t j = 0; j <= b->degree; ++j) {
            out->coefficients[i+j] = out->coefficients[i+j] ^ gf4_mul(a->coefficients[i], b->coefficients[j]);
        }
    }
    for (size_t i = 0; i < out->capacity; ++i) {
        if (0 != out->coefficients[i]) {
            out->degree = i;
        }
    }
}

// division and modulo
gf4_poly_t gf4_poly_div_x_to_deg(gf4_poly_t * poly, size_t deg) {
#ifdef SAFE
    assert(NULL != poly);
	assert(poly->degree <= deg);
#endif // SAFE
    gf4_poly_t out;
    out.capacity = poly->capacity;
    out.coefficients = calloc(out.capacity, sizeof(gf4_t));
    assert(NULL != out.coefficients);

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
#ifdef SAFE
    assert(NULL != out);
	assert(NULL != poly);
	assert(poly->degree <= deg);
#endif // SAFE
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
#ifdef SAFE
	assert(NULL != poly);
	assert(poly->degree <= deg);
#endif // SAFE
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
#ifdef SAFE
    assert(!gf4_poly_is_zero(b));
#endif // SAFE

    if (a->degree < b->degree) {
#ifdef SAFE
        gf4_poly_init_zero_byref(div, 0);
#else
        gf4_poly_init_zero_byref(div, a->capacity);
#endif // SAFE
        *rem = gf4_poly_clone(a);
        return;
    } else {
        gf4_poly_t q, r;
#ifdef SAFE
        gf4_poly_init_zero_byref(&q, 0);
#else
        gf4_poly_init_zero_byref(&q, a->capacity);
#endif
        r = gf4_poly_clone(a);

        gf4_t b_lead = b->coefficients[b->degree];

        while (!gf4_poly_is_zero(&r) && r.degree >= b->degree) {
            size_t deg = r.degree - b->degree;
            gf4_t tmp = gf4_div(r.coefficients[r.degree], b_lead);
            gf4_poly_add_ax_to_deg_inplace(&q, tmp, deg);

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
void invert_slow(gf4_poly_t * maybe_inverse, gf4_poly_t * poly, gf4_poly_t * modulus) {
#ifdef SAFE
    assert(!gf4_poly_is_zero(modulus));
    assert(NULL != maybe_inverse);
#endif
    if (gf4_poly_is_zero(poly)) {
        maybe_inverse = NULL;
        return;
    }
    gf4_poly_t a = gf4_poly_clone(poly);
    gf4_poly_t b = gf4_poly_clone(modulus);
#ifdef SAFE
    gf4_poly_t t = gf4_poly_init_zero(0);
    gf4_poly_t s = gf4_poly_init_zero(0);
    gf4_poly_t d = gf4_poly_init_zero(0);
    gf4_poly_t r = gf4_poly_init_zero(0);
#else
    gf4_poly_t t = gf4_poly_init_zero(poly->capacity);
    gf4_poly_t s = gf4_poly_init_zero(poly->capacity);
    gf4_poly_t d = gf4_poly_init_zero(poly->capacity);
    gf4_poly_t r = gf4_poly_init_zero(poly->capacity);
#endif

    while(!gf4_poly_is_zero(&b)) {
        gf4_poly_div_rem(&d, &r, &a, &b);

        // (a, b) <- (b, r)
        gf4_poly_copy_from(&a, &b); // a <- b
        gf4_poly_copy_from(&b, &r); // b <- r

        // (t, s) <- (s, t + d * s)
        gf4_poly_copy_from(&r, &t); // r <- t
        gf4_poly_copy_from(&t, &s); // t <- s
        gf4_poly_t tmp = gf4_poly_mul(&d, &s);
        tmp = gf4_poly_add(&r, &tmp);
        gf4_poly_copy_from(&s, &tmp); // s <- t + d * s

        for (size_t i = 0; i <= d.degree; ++i) {
            d.coefficients[i] = 0;
        }
        for (size_t i = 0; i <= r.degree; ++i) {
            r.coefficients[i] = 0;
        }
    }
    if(a.degree > 0) {
        maybe_inverse = NULL;
        gf4_poly_deinit(&t);
    } else {
        // t / a
        size_t deg = 0;
        for (size_t i = 0; i <= t.degree; ++i) {
            t.coefficients[i] = gf4_div(t.coefficients[i], a.coefficients[0]);
            if (0 != t.coefficients[i]) {
                deg = i;
            }
        }
        t.degree = deg;
        *maybe_inverse = t;
    }
    gf4_poly_deinit(&a);
    gf4_poly_deinit(&b);
    gf4_poly_deinit(&s);
    gf4_poly_deinit(&d);
    gf4_poly_deinit(&r);
}

// properties
bool gf4_poly_is_zero(gf4_poly_t * poly) {
    return 0 == poly->degree && 0 == poly->coefficients[0];
}



// helpers
void gf4_poly_pretty_print(gf4_poly_t * poly) {
    if (0 == poly->degree && 0 == poly->coefficients[0]) {
        printf("0\n");
    } else {
        char * sep = "";
        for (size_t i = 0; i <= poly->degree; ++i) {
            if (0 != poly->coefficients[i]) {
                printf("%s%s*x^%zu", sep, gf4_to_str(poly->coefficients[i]), i);
            }
        }
        printf("\n");
    }
}

gf4_poly_t gf4_poly_clone(gf4_poly_t * poly) {
    gf4_poly_t clone;
    clone.degree = poly->degree;
    clone.capacity = poly->capacity;
    clone.coefficients = calloc(clone.capacity, sizeof(gf4_t));
    assert(NULL != clone.coefficients);
    memcpy(clone.coefficients, poly->coefficients, clone.degree+1);
}