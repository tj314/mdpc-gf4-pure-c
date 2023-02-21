#ifndef MDPC_GF4_CONTEXTS_H
#define MDPC_GF4_CONTEXTS_H

#include <stdlib.h>
#include "gf4.h"
#include "gf4_poly.h"
#include "random.h"

typedef struct {
    gf4_poly_t second_block_G;
    size_t block_size;
} encoding_context_t;

typedef struct {
    gf4_poly_t h0;
    gf4_poly_t h1;
    size_t block_size;
} decoding_context_t;

/**
 * Generate contexts for encoding and decoding.
 *
 * generate polynomials h0 and h1 st. hamming weight of h0 (and also h1) is equal to block_weight.
 * generated polynomial h1 will be invertible mod (x^block_size + 1). h0 and h1 are used for decoding.
 * polynomial for encoding is calculated as follows: (h0 * h1^-1) mod (x^block_size + 1).
 *
 * This function will allocate memory for the polynomials. Do not initialize polynomials in out_enc_ctx and
 * out_dec_ctx yourself!
 *
 * @param out_enc_ctx memory location to store encoding polynomials and parameters to
 * @param out_dec_ctx memory location to store decoding polynomials and parameters to
 * @param block_size size of the circulant block of the matrices H, G
 * @param block_weight hamming weight of each row/columns of the circulant blocks of H and G
 */
void contexts_init(encoding_context_t * out_enc_ctx, decoding_context_t * out_dec_ctx, size_t block_size, size_t block_weight);

/**
 * Deinit contexts.
 *
 * Deinitializes polynomials contained within the contexts.
 *
 * @param enc_ctx memory location of the encoding context
 * @param dec_ctx memory location of the decoding context
 */
void contexts_deinit(encoding_context_t * enc_ctx, decoding_context_t * dec_ctx);
#endif //MDPC_GF4_CONTEXTS_H
