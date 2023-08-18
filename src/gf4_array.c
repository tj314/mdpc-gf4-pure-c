/**
 *  @file   gf4_array.c
 *  @brief  Arrays over GF(4).
 *  @author Tomáš Vavro
 *  @date   2023-08-12
 ***********************************************/

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

#include "gf4_array.h"

gf4_array_t gf4_array_init(size_t capacity, bool zero_out_new_memory) {
    assert(0 < capacity);
    gf4_array_t out;
    if (zero_out_new_memory) {
        out.array = calloc(capacity, sizeof(gf4_t));
    } else {
        out.array = malloc(capacity * sizeof(gf4_t));
    }
    if (NULL == out.array) {
        fprintf(stderr, "%s: Allocation error!\n", __func__);
        exit(-1);
    }
    out.capacity = capacity;
    return out;
}

gf4_array_t gf4_array_clone(gf4_array_t * in_array) {
    gf4_array_t clone;
    clone.capacity = in_array->capacity;
    clone.array = calloc(clone.capacity, sizeof(gf4_t));
    if (NULL == clone.array) {
        fprintf(stderr, "%s: Memory allocation failed!\n", __func__);
        exit(-1);
    }
    memcpy(clone.array, in_array->array, sizeof(gf4_t) * clone.capacity);
    return clone;
}

void gf4_array_resize(gf4_array_t * array, size_t new_capacity, bool zero_out_new_memory) {
    assert(NULL != array);
    assert(0 < new_capacity);
    gf4_t * tmp = realloc(array->array, new_capacity * sizeof(gf4_t));
    if (NULL == tmp) {
        free(array->array);
        fprintf(stderr, "gf4_array_resize: Memory reallocation failed!\n");
        exit(-1);
    }
    if (zero_out_new_memory) {
        for (size_t i = array->capacity; i < new_capacity; ++i) {
            tmp[i] = 0;
        }
    }
    array->array = tmp;
    array->capacity = new_capacity;
}

void gf4_array_deinit(gf4_array_t * array) {
    assert(NULL != array);
    free(array->array);
    array->array = NULL;
    array->capacity = 0;
}

// properties
size_t gf4_array_hamming_weight(gf4_array_t *array) {
    assert(NULL != array);
    size_t weight = 0;
    for (size_t i = 0; i < array->capacity; ++i) {
        weight += (0 != array->array[i]);
    }
    return weight;
}

// helpers
void gf4_array_print(gf4_array_t * array, FILE * stream, const char * end) {
    assert(NULL != array);
    assert(NULL != array->array);
    assert(NULL != stream);
    assert(NULL != end);
    for (size_t i = 0; i < array->capacity; ++i) {
        fprintf(stream, "%u ", array->array[i]);
    }
    fprintf(stream, "%s", end);
}

gf4_t gf4_array_sum(gf4_array_t * array) {
    assert(NULL != array);
    assert(NULL != array->array);
    gf4_t sum = 0;
    for (size_t i = 0; i < array->capacity; ++i) {
        sum = gf4_add(sum, array->array[i]);
    }
    return sum;
}

void gf4_array_zero_out(gf4_array_t * array) {
    assert(NULL != array);
    assert(NULL != array->array);
    for (size_t i = 0; i < array->capacity; ++i) {
        array->array[i] = 0;
    }
}