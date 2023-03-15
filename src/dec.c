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

bool dec_decode_symbol_flipping(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    gf4_poly_t syndrome_copy = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    size_t syndrome_weight = utils_hamming_weight(&syndrome);
    if (0 == syndrome_weight) {
        gf4_poly_deinit(&syndrome_copy);
        gf4_poly_deinit(&syndrome);
        return true;
    } else {
        for (size_t i = 0; i < num_iterations; ++i) {
            size_t min_weight = ctx->block_size + 1;
            size_t pos = 0;
            gf4_t a_max = 0;
            for (size_t j = 0; j < 2*ctx->block_size; ++j) {
                gf4_poly_copy(&syndrome_copy, &syndrome);
                for (gf4_t a = 1; a < GF4_MAX_VALUE; ++a) {
                    // s = s - ah_j, where h_j is j-th column of H
                    gf4_poly_t *h_block;
                    size_t actual_j;
                    if (j < ctx->block_size) {
                        h_block = &ctx->h0;
                        actual_j = j;
                    } else {
                        h_block = &ctx->h1;
                        actual_j = j - ctx->block_size;
                    }
                    size_t x = actual_j;
                    size_t idx = 0;
                    do {
                        syndrome_copy.coefficients[idx] ^= gf4_mul(h_block->coefficients[x], a);
                        x = (0 == x) ? ctx->block_size - 1: x-1;
                        ++idx;
                    } while (x != actual_j);

                    size_t new_syn_weight = utils_hamming_weight(&syndrome_copy);
                    if (new_syn_weight < min_weight) {
                        min_weight = new_syn_weight;
                        a_max = a;
                        pos = j;
                    }
                }
            }
            gf4_poly_t *h_block;
            size_t h_pos;
            if (pos < ctx->block_size) {
                h_block = &ctx->h0;
                h_pos = pos;
            } else {
                h_block = &ctx->h1;
                h_pos = pos - ctx->block_size;
            }
            size_t x = h_pos;
            size_t idx = 0;
            do {
                syndrome.coefficients[idx] ^= gf4_mul(h_block->coefficients[x], a_max);
                x = (0 == x) ? ctx->block_size - 1: x-1;
                ++idx;
            } while (x != h_pos);

            maybe_decoded->coefficients[pos] ^= a_max;
            if (0 == utils_hamming_weight(&syndrome)) {
                gf4_poly_deinit(&syndrome_copy);
                gf4_poly_deinit(&syndrome);
                return true;
            }
        }
        gf4_poly_deinit(&syndrome_copy);
        gf4_poly_deinit(&syndrome);
        return false;
    }
}

bool dec_decrypt(gf4_poly_t * out_decrypted, gf4_poly_t * in_encrypted, bool (*decode)(gf4_poly_t*, gf4_poly_t *, size_t, decoding_context_t *), size_t num_iterations, decoding_context_t * ctx) {
    bool decoding_success = decode(out_decrypted, in_encrypted, num_iterations, ctx);
    if (decoding_success) {
        memset(out_decrypted->coefficients + ctx->block_size, 0, ctx->block_size);
        return true;
    } else {
        gf4_poly_zero_out(out_decrypted);
        return false;
    }
}