//
// Created by tj on 10.8.2023.
//

#include "gf4_vector.h"

gf4_vector_t gf4_vector_init(size_t capacity) {
    assert(0 < capacity);
    gf4_vector_t out;
    out.capacity = capacity;
    out.length = 0;
    out.array = calloc(capacity, sizeof(gf4_t));
    if (NULL == out.array) {
        fprintf(stderr, "gf4_vector_init: Allocation error!\n");
        exit(-1);
    }
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