#include "enc.h"


void enc_encode(gf4_poly_t * out_encoded, gf4_poly_t * in_message, encoding_context_t * ctx) {
    assert(NULL != out_encoded);
    assert(NULL != in_message);
    assert(NULL != ctx);
    assert(out_encoded->capacity >= 2*ctx->block_size);
    assert(in_message >= ctx->block_size);
    memcpy(out_encoded->coefficients, in_message->coefficients, ctx->block_size);
    for (size_t i = 0; i < ctx->block_size; ++i) {
        gf4_t tmp = 0;
        for (size_t j = 0; j < ctx->block_size; ++j) {
            tmp ^= gf4_mul(in_message->coefficients[j], ctx->second_block_G.coefficients[(i + j) % ctx->block_size]);
        }
        out_encoded->coefficients[ctx->block_size + i] = tmp;
    }
}