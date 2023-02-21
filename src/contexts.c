#include "contexts.h"


// helper function
size_t poly_sum(gf4_poly_t * poly) {
    assert(NULL != poly);
    size_t sum = 0;
    for (size_t i = 0; i <= poly->degree; ++i) {
        sum += poly->coefficients[i];
    }
    return sum;
}

void generate_contexts(encoding_context_t * out_enc_ctx, decoding_context_t * out_dec_ctx, size_t block_size, size_t block_weight) {
    assert(NULL != out_enc_ctx);
    assert(NULL != out_dec_ctx);
    assert(block_weight <= block_size);
    size_t capacity = block_size + 1;
    gf4_poly_t modulus = gf4_poly_init_zero(capacity);
    gf4_poly_set_coefficient(&modulus, 0, 1);
    gf4_poly_set_coefficient(&modulus, block_size, 1);

    gf4_poly_t h0 = gf4_poly_init_zero(capacity);
    gf4_poly_t h1 = gf4_poly_init_zero(capacity);
    gf4_poly_t maybe_inverse = gf4_poly_init_zero(capacity);
    random_weighted_gf4_poly(&h0, block_weight);

    while (true) {
        while (0 == poly_sum(&h1)) {
            random_weighted_gf4_poly(&h1, block_weight);
        }
        bool inverted = gf4_poly_invert_slow(&maybe_inverse, &h1, &modulus);
        if (inverted) {
            gf4_poly_t tmp = gf4_poly_init_zero(capacity);
            gf4_poly_t div = gf4_poly_init_zero(capacity);
            gf4_poly_t rem = gf4_poly_init_zero(capacity);
            gf4_poly_mul(&tmp, &h1, &maybe_inverse);
            gf4_poly_div_rem(&div, &rem, &tmp, &modulus);
            bool correct_inverse = 0 == rem.degree && 1 == rem.coefficients[0];
            if (!correct_inverse) {
                // WTF???? this means invert function is incorrectly implemented
                gf4_poly_deinit(&tmp);
                gf4_poly_deinit(&rem);
                gf4_poly_deinit(&div);
                gf4_poly_deinit(&modulus);
                gf4_poly_deinit(&maybe_inverse);
                gf4_poly_deinit(&h0);
                gf4_poly_deinit(&h1);
                exit(-1);
            }

            // second_block_G_poly = (h0_poly * inverse) % modulus ;
            gf4_poly_zero_out(&tmp);
            gf4_poly_zero_out(&rem);
            gf4_poly_zero_out(&div);
            gf4_poly_mul(&tmp, &h0, &maybe_inverse);
            gf4_poly_div_rem(&div, &rem, &tmp, &modulus);
            gf4_poly_deinit(&tmp);
            gf4_poly_deinit(&rem);
            gf4_poly_deinit(&div);
            out_enc_ctx->block_size = block_size;
            out_enc_ctx->second_block_G = rem;
            out_dec_ctx->block_size = block_size;
            out_dec_ctx->h0 = h0;
            out_dec_ctx->h1 = h1;
            gf4_poly_deinit(&modulus);
            gf4_poly_deinit(&maybe_inverse);
            return;
        }
    }
}