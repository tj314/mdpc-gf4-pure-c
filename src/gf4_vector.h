/**
 *  @file   gf4_vector.h
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
gf4_vector_t gf4_vector_init(size_t capacity, bool zero_out_new_memory);

gf4_vector_t gf4_vector_init_with_length(size_t capacity, size_t length, bool zero_out_new_memory);

void gf4_vector_resize(gf4_vector_t * vector, size_t new_capacity, bool zero_out_new_memory);

void gf4_vector_deinit(gf4_vector_t * vector);

// properties
/**
 * @brief Find Hamming weight of a vector.
 *
 * vector must be initialized beforehand.
 *
 * @param vector pointer to a vector
 * @return Hamming weight of vector
 */
size_t gf4_vector_hamming_weight(gf4_vector_t *vector);

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
