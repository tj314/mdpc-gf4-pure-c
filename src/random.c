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

void random_gf4_vector(gf4_vector_t *vector, size_t size) {
    assert(NULL != vector);
    assert(size <= vector->capacity);
    for (size_t i = 0; i < size; ++i) {
        vector->array[i] = (gf4_t) random_from_range(0, GF4_MAX_VALUE);
    }
    vector->length = 1;
    for (size_t i = size - 1; i != 0; --i) {
        if (0 != vector->array[i]) {
            vector->length = i + 1;
            break;
        }
    }
}

void random_weighted_gf4_vector(gf4_vector_t *vector, size_t size, size_t weight) {
    assert(NULL != vector);
    assert(size <= vector->capacity);
    assert(weight <= size);
    for (size_t i = 0; i < weight; ++i) {
        vector->array[i] = (gf4_t) random_from_range(1, GF4_MAX_VALUE);
    }
    for (size_t i = 0; i < size - 1; ++i) {
        size_t j = (size_t)random_from_range(i, size - 1);
        gf4_t tmp = vector->array[i];
        vector->array[i] = vector->array[j];
        vector->array[j] = tmp;
    }
    vector->length = 1;
    for (size_t i = size - 1; i != 0; --i) {
        if (0 != vector->array[i]) {
            vector->length = i + 1;
            break;
        }
    }
}

void random_weighted_gf4_pairs_common(gf4_poly_t * poly, size_t size, size_t weight, size_t distance, gf4_t first, gf4_t second) {
    assert(NULL != poly);
    assert(poly->capacity >= size);
    assert(size > weight);
    size_t num_pairs = weight / 2;
    for (size_t i = 0; i < num_pairs; ++i) {
        size_t index = random_from_range(0, size - 1);
        size_t second_index = (index + distance) % size;
        while (poly->array[index] > 0 || poly->array[second_index] > 0) {
            index = random_from_range(0, size - 1);
            second_index = (index + distance) % size;
        }
        poly->array[index] = first;
        poly->array[second_index] = second;
    }
    gf4_poly_adjust_degree(poly, size - 1);
}

void random_weighted_gf4_vector_pairs_of_ones(gf4_vector_t *vector, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_pairs_common(vector, size, weight, distance, 1, 1);
}

void random_weighted_gf4_vector_pairs_of_one_alpha(gf4_vector_t *vector, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_pairs_common(vector, size, weight, distance, 1, 2);
}

void random_weighted_gf4_vector_pairs_of_alpha_one(gf4_vector_t *vector, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_pairs_common(vector, size, weight, distance, 2, 1);
}