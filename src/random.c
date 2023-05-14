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

#include "random.h"

void random_init() {
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        time_t t;
        srand((unsigned)time(&t));
    }
}

void random_force_reseed() {
    time_t t;
    srand((unsigned)time(&t));
}

size_t random_from_range(size_t low_bound_inclusive, size_t top_bound_inclusive) {
    assert(low_bound_inclusive < top_bound_inclusive);
    random_init();
    return low_bound_inclusive + (rand() % (top_bound_inclusive - low_bound_inclusive + 1));
}

void random_gf4_poly(gf4_poly_t * poly, size_t size) {
    assert(NULL != poly);
    assert(size <= poly->capacity);
    for (size_t i = 0; i < size; ++i) {
        poly->coefficients[i] = (gf4_t) random_from_range(0, GF4_MAX_VALUE);
    }
    poly->degree = 0;
    for (size_t i = size - 1; i != 0; --i) {
        if (0 != poly->coefficients[i]) {
            poly->degree = i;
            break;
        }
    }
}

void random_weighted_gf4_poly(gf4_poly_t * poly, size_t size, size_t weight) {
    assert(NULL != poly);
    assert(size <= poly->capacity);
    assert(weight <= size);
    for (size_t i = 0; i < weight; ++i) {
        poly->coefficients[i] = (gf4_t) random_from_range(1, GF4_MAX_VALUE);
    }
    for (size_t i = 0; i < size - 1; ++i) {
        size_t j = (size_t)random_from_range(i, size - 1);
        gf4_t tmp = poly->coefficients[i];
        poly->coefficients[i] = poly->coefficients[j];
        poly->coefficients[j] = tmp;
    }
    poly->degree = 0;
    for (size_t i = size - 1; i != 0; --i) {
        if (0 != poly->coefficients[i]) {
            poly->degree = i;
            break;
        }
    }
}

void random_weighted_gf4_pairs_common(gf4_poly_t * poly, size_t size, size_t weight, size_t distance, gf4_t first, gf4_t second) {
    assert(NULL != poly);
    size_t num_pairs = weight / 2;
    for (size_t i = 0; i < num_pairs; ++i) {
        size_t index = random_from_range(0, size - 1);
        while (poly->coefficients[index] > 0) {
            index = random_from_range(0, size - 1);
        }
        poly->coefficients[index] = first;
        poly->coefficients[(index + distance) % size] = second;
    }
    gf4_poly_adjust_degree(poly, size - 1);
}

void random_weighted_gf4_poly_pairs_of_ones(gf4_poly_t * poly, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_pairs_common(poly, size, weight, distance, 1, 1);
}

void random_weighted_gf4_poly_pairs_of_one_alpha(gf4_poly_t * poly, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_pairs_common(poly, size, weight, distance, 1, 2);
}