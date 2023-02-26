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

/**
 * Perform basic symbol-flipping decoding.
 *
 * maybe_decoded and in_vector must be initialized in advance and must have capacity = 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * @param maybe_decoded pointer to a polynomial tha will be used to store the decoded message
 * @param in_vector pointer to a polynomial representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx);


bool dec_decrypt(gf4_poly_t * out_decrypted, gf4_poly_t * in_encrypted, bool (*decode)(gf4_poly_t*, gf4_poly_t *, size_t, decoding_context_t *), size_t num_iterations, decoding_context_t * ctx);

#endif //MDPC_GF4_DEC_H
