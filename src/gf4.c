#include "gf4.h"

// helper function
bool gf4_is_in_range(gf4_t a) {
    return a <= GF4_MAX_VALUE;
}

const char * gf4_to_str(gf4_t a) {
    assert(gf4_is_in_range(a));
    switch (a) {
        case 0:
            return "0";
        case 1:
            return "1";
        case 2:
            return "a";
        case 3:
            return "(a+1)";
        default:
            return "";
    }
}


// operations
gf4_t gf4_add(gf4_t a, gf4_t b) {
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
    return a ^ b;
}

gf4_t gf4_mul(gf4_t a, gf4_t b) {
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
    return GF4_MULTIPLICATION[a][b];
}

gf4_t gf4_div(gf4_t a, gf4_t b) {
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
	assert(0 != b);
    return GF4_DIVISION[a][b-1];
}