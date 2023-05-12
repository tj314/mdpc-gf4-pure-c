/**
 *  @file   gf4.h
 *  @brief  Implementation of GF(4).
 *  @author Tomáš Vavro
 *  @date   2023-05-12
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

#ifndef GF4_GF4_H
#define GF4_GF4_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

#define GF4_MAX_VALUE 3 // inclusive!

typedef uint8_t gf4_t;


// cayley tables
static const gf4_t GF4_MULTIPLICATION[4][4] = {
        {0, 0, 0, 0},
        {0, 1, 2, 3},
        {0, 2, 3, 1},
        {0, 3, 1, 2}
};

static const gf4_t GF4_DIVISION[4][3] ={
        {0, 0, 0},
        {1, 3, 2},
        {2, 1, 3},
        {3, 2, 1}
};


// helper function

/**
 * @brief Check whether a is a valid value for a GF4 element.
 *
 * @param a possible GF4 value
 * @return true if 0 <= a <= GF4_MAX_VALUE, false otherwise
 */
bool gf4_is_in_range(gf4_t a);

/**
 * @brief Convert a GF4 element to string for printing.
 *
 * 0 -> "0"
 * 1 -> "1"
 * 2 -> "a"
 * 3 -> "(a+1)"
 *
 * @param a GF4 value to be converted
 * @return corresponding string representation
 */
const char * gf4_to_str(gf4_t a);

// operations

/**
 * @brief Compute a+b in GF4.
 *
 * In GF4, the following holds: a+b = a-b = a XOR b
 *
 * @param a GF4 element
 * @param b GF4 element
 * @return a+b
 */
gf4_t gf4_add(gf4_t a, gf4_t b);

/**
 * @brief Compute a*b in GF4.
 *
 * @param a GF4 element
 * @param b GF4 element
 * @return a*b
 */
gf4_t gf4_mul(gf4_t a, gf4_t b);

/**
 * @brief Compute a/b in GF4.
 *
 * Value b must be nonzero!
 *
 * @param a GF4 element
 * @param b nonzero GF4 element
 * @return a/b
 */
gf4_t gf4_div(gf4_t a, gf4_t b);

#endif // GF4_GF4_H
