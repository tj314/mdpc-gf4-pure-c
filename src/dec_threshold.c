/*
 This file is part of QC-MDPC McEliece over GF(4) implementation.
 Copyright (C) 2023 Tomáš Vavro

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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

bool dec_decode_symbol_flipping_threshold(gf4_array_t *maybe_decoded, gf4_array_t *in_array, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_array);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_array->capacity >= 2 * ctx->block_size);
    assert(NULL != ctx->threshold);

    gf4_array_t syndrome = gf4_array_init(ctx->block_size, true);
    dec_calculate_syndrome(&syndrome, in_array, ctx);
    memcpy(maybe_decoded->array, in_array->array, sizeof(gf4_t)*in_array->capacity);

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long) gf4_array_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            gf4_array_deinit(&syndrome);
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
                maybe_decoded->array[j] ^= a_max;
            }
        }
        gf4_array_zero_out(&syndrome);
        dec_calculate_syndrome(&syndrome, maybe_decoded, ctx);
    }
    gf4_array_deinit(&syndrome);
    ctx->elapsed_iterations = num_iterations;
    return false;
}