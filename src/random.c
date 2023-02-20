#include "random.h"

void random_init() {
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        time_t t;
        srand((unsigned)time(&t));
    }
}

size_t random_from_range(size_t low_bound_inclusive, size_t top_bound_inclusive) {
    assert(low_bound_inclusive < top_bound_inclusive);
    random_init();
    return low_bound_inclusive + (rand() % (top_bound_inclusive - low_bound_inclusive + 1));
}

void random_weighted_gf4_poly(gf4_poly_t * poly, size_t weight) {
    assert(NULL != poly);
    assert(weight <= poly->capacity);
    for (size_t i = 0; i < weight; ++i) {
        poly->coefficients[i] = (gf4_t) random_from_range(1, GF4_MAX_VALUE);
    }
    for (size_t i = 0; i < poly->capacity; ++i) {
        size_t j = (size_t)random_from_range(i, poly->capacity - 1);
        gf4_t tmp = poly->coefficients[i];
        poly->coefficients[i] = poly->coefficients[j];
        poly->coefficients[j] = tmp;
    }
    poly->degree = 0;
    for (size_t i = poly->capacity - 1; i != 0; --i) {
        if (0 != poly->coefficients[i]) {
            poly->degree = i;
            break;
        }
    }
}