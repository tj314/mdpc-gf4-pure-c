//
// Created by tj on 10.8.2023.
//

#ifndef MDPC_GF4_GF4_VECTOR_H
#define MDPC_GF4_GF4_VECTOR_H

#include <stdlib.h>
#include <stdbool.h>
#include "gf4.h"
/**
 * @brief Structure that represents a vector over GF4.
 *
 * It holds that capacity >= length.
 */
typedef struct {
    gf4_t * array; ///< an array
    size_t length; ///< actual used length of the array
    size_t capacity; ///< allocated amount of memory
} gf4_vector_t;

// initialization
gf4_vector_t gf4_vector_init(size_t capacity);

gf4_vector_t gf4_vector_init_with_length(size_t capacity, size_t length);

void gf4_vector_resize(gf4_vector_t * vector, size_t new_capacity, bool zero_out_new_memory);

void gf4_vector_deinit(gf4_vector_t * vector);

// helpers
/**
 * @brief Print array up to max.
 *
 * vector must be initialized beforehand.
 * max must be less or equal to vector->capacity
 * e.g.: [0, 1, 0, 2, 2, 1, 0], max=4 --> [0, 1, 0, (a+1)]
 *
 * @param vector pointer to a vector to be printed
 * @param max integer, amount of items to be printed
 * @param stream stream to be used (e.g. stdout, stderr...)
 * @param end last character to be printed (e. g. line ending or space)
 */
void gf4_vector_print(gf4_vector_t * vector, size_t max, FILE * stream, const char * end);

#endif //MDPC_GF4_GF4_VECTOR_H
