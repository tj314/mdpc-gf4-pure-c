#include "dec.h"

void dec_calculate_syndrome(gf4_poly_t * out_syndrome, gf4_poly_t * in_vector, decoding_context_t * ctx) {
    assert(NULL != out_syndrome);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(out_syndrome->capacity >= ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);
    for (size_t i = ctx->block_size; i > 0; --i) {
        gf4_t tmp = 0;
        for (size_t j = 0; j < ctx->block_size; ++j) {
            tmp ^= gf4_mul(ctx->h0.coefficients[(i+j) % ctx->block_size], in_vector->coefficients[j]);
            tmp ^= gf4_mul(ctx->h1.coefficients[(i+j) % ctx->block_size], in_vector->coefficients[ctx->block_size + j]);
        }
        out_syndrome->coefficients[ctx->block_size - i] = tmp;
    }
    out_syndrome->degree = 0;
    for (size_t i = out_syndrome->capacity; i > 0; --i) {
        if (0 != out_syndrome->coefficients[i]) {
            out_syndrome->degree = i;
            break;
        }
    }
}