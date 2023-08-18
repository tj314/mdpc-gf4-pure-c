/**
 *  @file   gf4_array.h
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

#ifndef MDPC_GF4_GF4_ARRAY_H
#define MDPC_GF4_GF4_ARRAY_H

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "gf4.h"

/**
 * @brief Structure that represents a array over GF4.
 */
typedef struct {
    gf4_t * array; ///< an array
    size_t capacity; ///< allocated amount of memory
} gf4_array_t;

// initialization
gf4_array_t gf4_array_init(size_t capacity, bool zero_out_new_memory);

gf4_array_t gf4_array_clone(gf4_array_t * in_array);

void gf4_array_resize(gf4_array_t * array, size_t new_capacity, bool zero_out_new_memory);

void gf4_array_deinit(gf4_array_t * array);

// properties
/**
 * @brief Find Hamming weight of an array.
 *
 * array must be initialized beforehand.
 *
 * @param array pointer to an array
 * @return Hamming weight of array
 */
size_t gf4_array_hamming_weight(gf4_array_t *array);

// helpers
/**
 * @brief Print array up to max.
 *
 * array must be initialized beforehand.
 * max must be less or equal to array->capacity
 * e.g.: [0, 1, 0, 2, 2, 1, 0], max=4 --> [0, 1, 0, (a+1)]
 *
 * @param array pointer to an array to be printed
 * @param stream stream to be used (e.g. stdout, stderr...)
 * @param end last character to be printed (e. g. line ending or space)
 */
void gf4_array_print(gf4_array_t * array, FILE * stream, const char * end);

/**
 * @brief Calculate the sum of all the elements in the array.
 *
 * @param array a pointer to an array to be summed
 * @return sum of all elements in the array
 */
gf4_t gf4_array_sum(gf4_array_t * array);

void gf4_array_zero_out(gf4_array_t * array);

#endif //MDPC_GF4_GF4_ARRAY_H
