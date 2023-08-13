/**
 *  @file   gf4_vector.c
 *  @brief  Vectors over GF(4).
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

#include "gf4_vector.h"

gf4_vector_t gf4_vector_init(size_t capacity, bool zero_out_new_memory) {
    assert(0 < capacity);
    gf4_vector_t out;
    if (zero_out_new_memory) {
        out.array = calloc(capacity, sizeof(gf4_t));
    } else {
        out.array = malloc(capacity * sizeof(gf4_t));
    }
    if (NULL == out.array) {
        fprintf(stderr, "gf4_vector_init: Allocation error!\n");
        exit(-1);
    }
    out.capacity = capacity;
    return out;
}

gf4_vector_t gf4_vector_init_with_length(size_t capacity, size_t length, bool zero_out_new_memory) {
    assert(length <= capacity);
    gf4_vector_t out = gf4_vector_init(capacity, zero_out_new_memory);
    out.length = length;
    return out;
}

void gf4_vector_resize(gf4_vector_t * vector, size_t new_capacity, bool zero_out_new_memory) {
    assert(NULL != vector);
    assert(vector->length <= new_capacity);
    gf4_t * tmp = realloc(vector->array, new_capacity*sizeof(gf4_t));
    if (NULL == tmp) {
        free(vector->array);
        fprintf(stderr, "gf4_vector_resize: Memory reallocation failed!\n");
        exit(-1);
    }
    if (zero_out_new_memory) {
        for (size_t i = vector->capacity; i < new_capacity; ++i) {
            tmp[i] = 0;
        }
    }
    vector->array = tmp;
    vector->capacity = new_capacity;
}

void gf4_vector_deinit(gf4_vector_t * vector) {
    assert(NULL != vector);
    free(vector->array);
    vector->array = NULL;
    vector->capacity = 0;
    vector->length = 0;
}

// properties
size_t gf4_vector_hamming_weight(gf4_vector_t *vector) {
    size_t weight = 0;
    for (size_t i = 0; i < vector->capacity; ++i) {
        weight += (0 != vector->array[i]);
    }
    return weight;
}

// helpers
void gf4_vector_print(gf4_vector_t * vector, size_t max, FILE * stream, const char * end) {
    assert(max <= vector->capacity);
    for (size_t i = 0; i < max; ++i) {
        fprintf(stream, "%u ", vector->array[i]);
    }
    fprintf(stream, "%s", end);
}