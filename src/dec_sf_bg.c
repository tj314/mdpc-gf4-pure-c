#include "dec.h"

void dec_decode_symbol_flipping_bg_iteration(gf4_poly_t * maybe_decoded, gf4_poly_t * syndrome, decoding_context_t * ctx) {
    long syndrome_weight = (long)utils_hamming_weight(syndrome);
    if (0 == syndrome_weight) {
        return;
    }
    long threshold = ctx->threshold(syndrome_weight);
    gf4_t * black = calloc(2*ctx->block_size, sizeof(gf4_t));
    assert(NULL != black);
    gf4_t * gray = calloc(2*ctx->block_size, sizeof(gf4_t));
    assert(NULL != gray);

    // iteration with setting black and gray
    size_t total_flipped = 0;
    for (size_t j = 0; j < 2*ctx->block_size; ++j) {
        gf4_poly_t *h_block;
        size_t actual_j;
        if (j < ctx->block_size) {
            h_block = &ctx->h0;
            actual_j = j;
        } else {
            h_block = &ctx->h1;
            actual_j = j - ctx->block_size;
        }
        long sigma_max = -1;
        gf4_t a_max = 0;
        for (gf4_t a = 1; a <= GF4_MAX_VALUE; ++a) {
            long sigma = dec_calculate_new_sigma(h_block, syndrome, syndrome_weight, a, actual_j, ctx);
            if (sigma > sigma_max) {
                sigma_max = sigma;
                a_max = a;
            }
        }

        if (sigma_max > threshold) {
            maybe_decoded->coefficients[j] ^= a_max;
            black[j] = a_max;
            total_flipped++;
        }
        long bound = (threshold < ctx->delta_setting) ? threshold - ctx->delta_setting : 0;
        if (sigma_max > bound) {
            gray[j] = a_max;
        }
    }
    //fprintf(stderr, "total flipped: %zu\n", total_flipped);
    gf4_poly_zero_out(syndrome);
    dec_calculate_syndrome(syndrome, maybe_decoded, ctx);
    syndrome_weight = (long) utils_hamming_weight(syndrome);

    // black
    // threshold = (long)(37.0/3);
    // threshold = (long)(((float)ctx->threshold(syndrome_weight)) / 2);
    threshold = 0;
    size_t black_flipped = 0;
    for (size_t j = 0; j < 2*ctx->block_size; ++j) {
        gf4_poly_t *h_block;
        size_t actual_j;
        if (j < ctx->block_size) {
            h_block = &ctx->h0;
            actual_j = j;
        } else {
            h_block = &ctx->h1;
            actual_j = j - ctx->block_size;
        }

        long sigma = dec_calculate_new_sigma(h_block, syndrome, syndrome_weight, black[j], actual_j, ctx);
        if (sigma > threshold) {
            maybe_decoded->coefficients[j] ^= black[j];
            black_flipped++;
        }
    }
    gf4_poly_zero_out(syndrome);
    dec_calculate_syndrome(syndrome, maybe_decoded, ctx);
    syndrome_weight = (long) utils_hamming_weight(syndrome);

    //fprintf(stderr, "black flipped: %zu\n", black_flipped);

    // gray
    // threshold = ctx->threshold(syndrome_weight);
    threshold = 9;
    size_t gray_flipped = 0;
    for (size_t j = 0; j < 2*ctx->block_size; ++j) {
        gf4_poly_t *h_block;
        size_t actual_j;
        if (j < ctx->block_size) {
            h_block = &ctx->h0;
            actual_j = j;
        } else {
            h_block = &ctx->h1;
            actual_j = j - ctx->block_size;
        }

        long sigma = dec_calculate_new_sigma(h_block, syndrome, syndrome_weight, gray[j], actual_j, ctx);
        if (sigma > threshold) {
            maybe_decoded->coefficients[j] ^= gray[j];
            gray_flipped++;
        }
    }
    gf4_poly_zero_out(syndrome);
    dec_calculate_syndrome(syndrome, maybe_decoded, ctx);
    //fprintf(stderr, "gray flipped: %zu\n", gray_flipped);

    free(black);
    free(gray);
}

bool dec_decode_symbol_flipping_bg(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);
    assert(NULL != ctx->threshold);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    if (gf4_poly_is_zero(&syndrome)) {
        gf4_poly_deinit(&syndrome);
        ctx->elapsed_iterations = 0;
        return true;
    }

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long)utils_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            gf4_poly_deinit(&syndrome);
            ctx->elapsed_iterations = i;
            return true;
        }
        dec_decode_symbol_flipping_bg_iteration(maybe_decoded, &syndrome, ctx);
    }
    gf4_poly_deinit(&syndrome);
    ctx->elapsed_iterations = num_iterations;
    return false;
}