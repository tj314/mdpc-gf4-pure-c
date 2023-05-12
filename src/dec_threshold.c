#include "dec.h"

long dec_calculate_threshold_0(long syndrome_weight) {
    // 0.0248577875*syndrome_weight - 29.1143817
    double tmp = 0.0248577875*(double)syndrome_weight - 29.1143817;
    if (tmp <= 0.0) tmp = 0.0;
    return (long)tmp;
}

long dec_calculate_threshold_1(long syndrome_weight) {
    // 0.0248577875*syndrome_weight - 29.1143817
    double tmp = 0.0248577875*(double)syndrome_weight - 28.1143817;
    if (tmp <= 0.0) tmp = 0.0;
    return (long)tmp;
}

long dec_calculate_threshold_2(long syndrome_weight) {
    // 0.0248577875*syndrome_weight - 29.1143817
    double tmp = 0.0248577875*(double)syndrome_weight - 27.1143817;
    if (tmp <= 0.0) tmp = 0.0;
    return (long)tmp;
}

long dec_calculate_threshold_3(long syndrome_weight) {
    // 0.0248577875*syndrome_weight - 29.1143817
    double tmp = 0.0248577875*(double)syndrome_weight - 26.1143817;
    if (tmp <= 0.0) tmp = 0.0;
    return (long)tmp;
}

long dec_calculate_threshold_4(long syndrome_weight) {
    // 0.0248577875*syndrome_weight - 29.1143817
    double tmp = 0.0248577875*(double)syndrome_weight - 25.1143817;
    if (tmp <= 0.0) tmp = 0.0;
    return (long)tmp;
}

long dec_calculate_threshold_5(long syndrome_weight) {
    // 0.0248577875*syndrome_weight - 29.1143817
    double tmp = 0.0248577875*(double)syndrome_weight - 24.1143817;
    if (tmp <= 0.0) tmp = 0.0;
    return (long)tmp;
}


bool dec_decode_symbol_flipping_threshold_recalculate_inplace(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);
    assert(NULL != ctx->threshold);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    long min_syndrome = -1;
    long max_syndrome = -1;
    size_t max_iter = 0;
    size_t min_iter = 0;

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long)utils_hamming_weight(&syndrome);
        if (min_syndrome < 0) {
            min_syndrome = syndrome_weight;
            max_syndrome = syndrome_weight;
            min_iter = i;
            max_iter = i;
        }
        if (0 == syndrome_weight) {
            gf4_poly_deinit(&syndrome);
            ctx->elapsed_iterations = i;
            fprintf(stderr, "elapsed iterations: %zu\n", i);
            return true;
        }
        long threshold = ctx->threshold(syndrome_weight);
        // fprintf(stderr, "syndrome weight = %ld, threshold = %ld\n", syndrome_weight, threshold);
        // fprintf(stderr, "Flipped positions: \n");
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
                long sigma = dec_calculate_new_sigma(h_block, &syndrome, syndrome_weight, a, actual_j, ctx);
                if (sigma > sigma_max) {
                    sigma_max = sigma;
                    a_max = a;
                }
            }

            if (sigma_max > threshold) {
                dec_flip_symbol(h_block, &syndrome, maybe_decoded, a_max, actual_j, j, ctx);
                syndrome_weight = (long)utils_hamming_weight(&syndrome);
                threshold = ctx->threshold(syndrome_weight);
                if (syndrome_weight < min_syndrome) {
                    min_syndrome = syndrome_weight;
                    min_iter = i;
                }
                if (syndrome_weight > max_syndrome) {
                    max_syndrome = syndrome_weight;
                    max_iter = i;
                }
            }
        }
    }
    gf4_poly_deinit(&syndrome);
    fprintf(stderr, "min=%ld at %zu max=%ld at %zu\n", min_syndrome, min_iter, max_syndrome, max_iter);

    ctx->elapsed_iterations = num_iterations;
    return false;
}

bool dec_decode_symbol_flipping_threshold(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);
    assert(NULL != ctx->threshold);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long)utils_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            gf4_poly_deinit(&syndrome);
            ctx->elapsed_iterations = i;
            return true;
        }
        long threshold = ctx->threshold(syndrome_weight);
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
                long sigma = dec_calculate_new_sigma(h_block, &syndrome, syndrome_weight, a, actual_j, ctx);
                if (sigma > sigma_max) {
                    sigma_max = sigma;
                    a_max = a;
                }
            }

            if (sigma_max > threshold) {
                maybe_decoded->coefficients[j] ^= a_max;
            }
        }
        gf4_poly_zero_out(&syndrome);
        dec_calculate_syndrome(&syndrome, maybe_decoded, ctx);
    }
    gf4_poly_deinit(&syndrome);
    ctx->elapsed_iterations = num_iterations;
    return false;
}