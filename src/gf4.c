#include "gf4.h"

// helper function
bool gf4_is_in_range(gf4_t a) {
    return a <= 3;
}

const char * gf4_to_str(gf4_t a) {
#ifdef SAFE
    assert(gf4_is_in_range(a));
#endif // SAFE
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
#ifdef SAFE
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
#endif // SAFE
    return a ^ b;
}

gf4_t gf4_mul(gf4_t a, gf4_t b) {
#ifdef SAFE
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
#endif // SAFE
    return GF4_MULTIPLICATION[a][b];
}

gf4_t gf4_div(gf4_t a, gf4_t b) {
#ifdef SAFE
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
	assert(0 != b);
#endif // SAFE
    return GF4_DIVISION[a][b-1];
}