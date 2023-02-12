#ifndef GF4_GF4_H 
#define GF4_GF4_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

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
bool gf4_is_in_range(gf4_t a);
const char * gf4_to_str(gf4_t a);


// operations
gf4_t gf4_add(gf4_t a, gf4_t b);
gf4_t gf4_mul(gf4_t a, gf4_t b);
gf4_t gf4_div(gf4_t a, gf4_t b);

#endif // GF4_GF4_H
