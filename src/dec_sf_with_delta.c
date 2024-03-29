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

bool dec_decode_symbol_flipping_delta(gf4_array_t *maybe_decoded, gf4_array_t *in_array, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_array);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_array->capacity >= 2 * ctx->block_size);
    assert(ctx->delta_setting >= 0);

    gf4_array_t syndrome = gf4_array_init(ctx->block_size, true);
    dec_calculate_syndrome(&syndrome, in_array, ctx);
    memcpy(maybe_decoded->array, in_array->array, sizeof(gf4_t)*in_array->capacity);
    ctx->elapsed_iterations = 0;

    long * sigmas = calloc(2 * ctx->block_size, sizeof(long));
    assert(NULL != sigmas);
    gf4_t * values = calloc(2 * ctx->block_size, sizeof(gf4_t));
    assert(NULL != values);

    const long DELTA = ctx->delta_setting;

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long) gf4_array_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            free(sigmas);
            free(values);
            gf4_array_deinit(&syndrome);
            ctx->elapsed_iterations = i;
            return true;
        }
        long sigma_max = -1;
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            sigmas[j] = -1;
            values[j] = 0;
            gf4_poly_t *h_block;
            size_t actual_j;
            if (j < ctx->block_size) {
                h_block = &ctx->h0;
                actual_j = j;
            } else {
                h_block = &ctx->h1;
                actual_j = j - ctx->block_size;
            }
            for (gf4_t a = 1; a <= GF4_MAX_VALUE; ++a) {
                // s = s - a*h_j, where h_j is j-th column of H
                size_t x = actual_j;
                size_t idx = 0;
                long new_syndrome_weight = 0; // hamming weight of s-a*h_j
                do {
                    gf4_t tmp = syndrome.array[idx] ^ gf4_mul(h_block->coefficients.array[x], a);
                    if (0 != tmp) {
                        new_syndrome_weight += 1;
                    }
                    x = (0 == x) ? ctx->block_size - 1 : x-1;
                    ++idx;
                } while (x != actual_j);

                long sigma = syndrome_weight - new_syndrome_weight;
                if (sigma > sigma_max) {
                    sigma_max = sigma;
                }
                if (sigma > sigmas[j]) {
                    sigmas[j] = sigma;
                    values[j] = a;
                }
            }
        }

        long bound = ((sigma_max - DELTA) >= 0) ? sigma_max - DELTA : 0;
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            if (sigmas[j] < bound) continue;
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
                syndrome.array[idx] ^= gf4_mul(h_block->coefficients.array[x], values[j]);
                x = (0 == x) ? ctx->block_size - 1 : x-1;
                ++idx;
            } while (x != actual_j);
            maybe_decoded->array[j] ^= values[j];
        }
    }
    free(sigmas);
    free(values);
    gf4_array_deinit(&syndrome);
    ctx->elapsed_iterations = num_iterations;
    return false;
}