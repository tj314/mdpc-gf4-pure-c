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

void random_gf4_array(gf4_array_t *array, size_t size) {
    assert(NULL != array);
    assert(size <= array->capacity);
    for (size_t i = 0; i < size; ++i) {
        array->array[i] = (gf4_t) random_from_range(0, GF4_MAX_VALUE);
    }
}

void random_weighted_gf4_array(gf4_array_t *array, size_t size, size_t weight) {
    assert(NULL != array);
    assert(size <= array->capacity);
    assert(weight <= size);
    // place weight nonzero entries on the first weight positions
    for (size_t i = 0; i < weight; ++i) {
        array->array[i] = (gf4_t) random_from_range(1, GF4_MAX_VALUE);
    }
    // shuffle
    for (size_t i = 0; i < size - 1; ++i) {
        size_t j = (size_t)random_from_range(i, size - 1);
        gf4_t tmp = array->array[i];
        array->array[i] = array->array[j];
        array->array[j] = tmp;
    }
}

void random_weighted_gf4_vector_pairs_common(gf4_array_t *vector, size_t size, size_t weight, size_t distance, gf4_t first, gf4_t second) {
    assert(NULL != vector);
    assert(vector->capacity >= size);
    assert(size > weight);
    size_t num_pairs = weight / 2;
    for (size_t i = 0; i < num_pairs; ++i) {
        size_t index = random_from_range(0, size - 1);
        size_t second_index = (index + distance) % size;
        while (vector->array[index] > 0 || vector->array[second_index] > 0) {
            index = random_from_range(0, size - 1);
            second_index = (index + distance) % size;
        }
        vector->array[index] = first;
        vector->array[second_index] = second;
    }
}

void random_weighted_gf4_array_pairs_of_ones(gf4_array_t *array, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_vector_pairs_common(array, size, weight, distance, 1, 1);
}

void random_weighted_gf4_array_pairs_of_one_alpha(gf4_array_t *array, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_vector_pairs_common(array, size, weight, distance, 1, 2);
}

void random_weighted_gf4_array_pairs_of_alpha_one(gf4_array_t *array, size_t size, size_t weight, size_t distance) {
    random_weighted_gf4_vector_pairs_common(array, size, weight, distance, 2, 1);
}