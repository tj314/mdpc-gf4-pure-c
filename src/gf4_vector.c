//
// Created by tj on 10.8.2023.
//

#include "gf4_vector.h"

gf4_vector_t gf4_vector_init(size_t capacity) {
    // no need to assert anything, gf4_vector_init_with_length will take care of it
    return gf4_vector_init_with_length(capacity, 0);
}

gf4_vector_t gf4_vector_init_with_length(size_t capacity, size_t length) {
    assert(0 < capacity);
    assert(length <= capacity);
    gf4_vector_t out;
    out.array = calloc(capacity, sizeof(gf4_t));
    if (NULL == out.array) {
        fprintf(stderr, "gf4_vector_init: Allocation error!\n");
        exit(-1);
    }
    out.capacity = capacity;
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