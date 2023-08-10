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

gf4_vector_t gf4_vector_init(size_t capacity);

void gf4_vector_resize(gf4_vector_t * vector, size_t new_capacity, bool zero_out_new_memory);

void gf4_vector_deinit(gf4_vector_t * vector);

#endif //MDPC_GF4_GF4_VECTOR_H
