#ifndef MDPC_GF4_DEC_H
#define MDPC_GF4_DEC_H

#include "gf4.h"
#include "gf4_poly.h"
#include "contexts.h"

/**
 * Calculate the syndrome of a vector.
 *
 * out_syndrome must be initialized in advance. It must have capacity of at least block_size. It must be a zero polynomial.
 *
 * @param out_syndrome pointer to a polynomial that will be used to store the calculated syndrome
 * @param in_vector pointer to a polynomial representing a vector
 * @param ctx a valid decoding context
 */
void dec_calculate_syndrome(gf4_poly_t * out_syndrome, gf4_poly_t * in_vector, decoding_context_t * ctx);

#endif //MDPC_GF4_DEC_H
