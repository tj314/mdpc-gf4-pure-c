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
 * Check whether a is a valid value for a GF4 element.
 *
 * @param a possible GF4 value
 * @return true if 0 <= a <= GF4_MAX_VALUE, false otherwise
 */
bool gf4_is_in_range(gf4_t a);

/**
 * Convert a GF4 element to string for printing.
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
 * Compute a+b in GF4.
 *
 * In GF4, the following holds: a+b = a-b = a XOR b
 *
 * @param a GF4 element
 * @param b GF4 element
 * @return a+b
 */
gf4_t gf4_add(gf4_t a, gf4_t b);

/**
 * Compute a*b in GF4.
 *
 * @param a GF4 element
 * @param b GF4 element
 * @return a*b
 */
gf4_t gf4_mul(gf4_t a, gf4_t b);

/**
 * Compute a/b in GF4.
 *
 * Value b must be nonzero!
 *
 * @param a GF4 element
 * @param b nonzero GF4 element
 * @return a/b
 */
gf4_t gf4_div(gf4_t a, gf4_t b);

#endif // GF4_GF4_H
