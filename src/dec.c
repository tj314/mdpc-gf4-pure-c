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
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    long syndrome_weight = (long)utils_hamming_weight(&syndrome);
    if (0 == syndrome_weight) {
        gf4_poly_deinit(&syndrome);
        return true;
    } else {
        for (size_t i = 0; i < num_iterations; ++i) {
            long sigma_max = -1;
            size_t pos = 0;
            gf4_t a_max = 0;
            for (size_t j = 0; j < 2*ctx->block_size; ++j) {
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
                    long tmp_hamming_weight = 0; // hamming weight of s-a*h_j
                    do {
                        gf4_t tmp = syndrome.coefficients[idx] ^ gf4_mul(h_block->coefficients[x], a);
                        if (0 != tmp) {
                            tmp_hamming_weight += 1;
                        }
                        x = (0 == x) ? ctx->block_size - 1: x-1;
                        ++idx;
                    } while (x != actual_j);

                    long sigma = syndrome_weight - tmp_hamming_weight;
                    if (sigma > sigma_max) {
                        sigma_max = sigma;
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
            // syndrome_weight = (long)utils_hamming_weight(&syndrome);
            if (0 == utils_hamming_weight(&syndrome)) {
                gf4_poly_deinit(&syndrome);
                return true;
            }
        }
        gf4_poly_deinit(&syndrome);
        return false;
    }
}

bool dec_decode_symbol_flipping_2(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    long syndrome_weight = (long)utils_hamming_weight(&syndrome);
    if (0 == syndrome_weight) {
        gf4_poly_deinit(&syndrome);
        return true;
    } else {
        for (size_t i = 0; i < num_iterations; ++i) {
            long max_new_syndrome_weight = (long)ctx->block_size + 1; // it is impossible for w(syndrome) to be more than block size
            size_t pos = 0;
            gf4_t a_max = 0;
            for (size_t j = 0; j < 2*ctx->block_size; ++j) {
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
                    long tmp_hamming_weight = 0; // hamming weight of s-a*h_j
                    do {
                        gf4_t tmp = syndrome.coefficients[idx] ^ gf4_mul(h_block->coefficients[x], a);
                        if (0 != tmp) {
                            tmp_hamming_weight += 1;
                        }
                        x = (0 == x) ? ctx->block_size - 1: x-1;
                        ++idx;
                    } while (x != actual_j);

                    if (tmp_hamming_weight < max_new_syndrome_weight) {
                        max_new_syndrome_weight = tmp_hamming_weight;
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
            // syndrome_weight = (long)utils_hamming_weight(&syndrome);
            if (0 == utils_hamming_weight(&syndrome)) {
                gf4_poly_deinit(&syndrome);
                return true;
            }
        }
        gf4_poly_deinit(&syndrome);
        return false;
    }
}

bool dec_decode_symbol_flipping_delta(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    // const long DELTA = 3; // TODO: use the best setting
    const long delta = ctx->delta_setting;

    long syndrome_weight = (long)utils_hamming_weight(&syndrome);
    if (0 == syndrome_weight) {
        gf4_poly_deinit(&syndrome);
        return true;
    } else {
        size_t num_positions_to_flip = 0;
        long * sigmas = malloc(2*ctx->block_size * sizeof(long));
        gf4_t * symbols_to_flip_to = malloc(2*ctx->block_size * sizeof(gf4_t));
        if (NULL == sigmas || NULL == symbols_to_flip_to) {
            fprintf(stderr, "dec_decode_symbol_flipping_delta: Allocation error!");
            exit(-1);
        }
        for (size_t i = 0; i < num_iterations; ++i) {
            long sigma_max = -1;
            memset(sigmas, -1, 2*ctx->block_size);
            memset(symbols_to_flip_to, 0, 2*ctx->block_size);
            for (size_t j = 0; j < 2*ctx->block_size; ++j) {
                for (gf4_t a = 1; a < GF4_MAX_VALUE; ++a) {
                    // s = s - ah_j, where h_j is j-th column of H
                    // first, find the correct block
                    gf4_poly_t *h_block;
                    size_t actual_j;
                    if (j < ctx->block_size) {
                        h_block = &ctx->h0;
                        actual_j = j;
                    } else {
                        h_block = &ctx->h1;
                        actual_j = j - ctx->block_size;
                    }
                    // second, calculate the hamming weight of s-a*h_j
                    // the h_j column is constructed by some clever indexing
                    size_t x = actual_j;
                    size_t idx = 0;
                    long tmp_hamming_weight = 0; // hamming weight of s-a*h_j
                    do {
                        gf4_t tmp = syndrome.coefficients[idx] ^ gf4_mul(h_block->coefficients[x], a);
                        if (0 != tmp) {
                            tmp_hamming_weight += 1;
                        }
                        x = (0 == x) ? ctx->block_size - 1: x-1;
                        ++idx;
                    } while (x != actual_j);

                    long sigma = syndrome_weight - tmp_hamming_weight;
                    if (sigma > sigma_max) {
                        sigma_max = sigma;
                    }
                    if (sigma > sigmas[j]) {
                        sigmas[j] = sigma;
                        symbols_to_flip_to[j] = a;
                    }
                }
            }
            // update the syndrome and maybe_decoded
            for (size_t j = 0; j < 2*ctx->block_size; ++j) {
                if (sigmas[j] > sigma_max - DELTA) {
                    gf4_poly_t *h_block;
                    size_t h_pos;
                    if (j < ctx->block_size) {
                        h_block = &ctx->h0;
                        h_pos = j;
                    } else {
                        h_block = &ctx->h1;
                        h_pos = j - ctx->block_size;
                    }
                    size_t x = h_pos;
                    size_t idx = 0;
                    do {
                        syndrome.coefficients[idx] ^= gf4_mul(h_block->coefficients[x], symbols_to_flip_to[j]);
                        x = (0 == x) ? ctx->block_size - 1: x-1;
                        ++idx;
                    } while (x != h_pos);

                    maybe_decoded->coefficients[j] ^= symbols_to_flip_to[j];
                }
            }
            // syndrome_weight = (long)utils_hamming_weight(&syndrome);
            if (0 == utils_hamming_weight(&syndrome)) {
                gf4_poly_deinit(&syndrome);
                return true;
            }
        }
        free(sigmas);
        free(symbols_to_flip_to);
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