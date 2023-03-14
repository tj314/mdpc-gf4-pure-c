#include "utils.h"

size_t utils_hamming_weight(gf4_poly_t * vector) {
    size_t weight = 0;
    for (size_t i = 0; i < vector->capacity; ++i) {
        weight += (0 != vector->coefficients[i]);
    }
    return weight;
}