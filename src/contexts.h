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


void contexts_genarate(encoding_context_t * out_enc_ctx, decoding_context_t * out_dec_ctx, size_t block_size, size_t block_weight);
#endif //MDPC_GF4_CONTEXTS_H
