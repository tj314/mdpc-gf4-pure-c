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

#ifndef WRITE_WEIGHTS
bool dec_decode_symbol_flipping(gf4_array_t *maybe_decoded, gf4_array_t *in_array, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_array);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_array->capacity >= 2 * ctx->block_size);

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
        long sigma_max = -1;
        size_t pos = 0;
        gf4_t a_max = 0;
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
            for (gf4_t a = 1; a <= GF4_MAX_VALUE; ++a) {
                // s = s - a*h_j, where h_j is j-th column of H
                size_t x = actual_j;
                size_t idx = 0;
                long tmp_hamming_weight = 0; // hamming weight of s-a*h_j
                do {
                    gf4_t tmp = syndrome.array[idx] ^ gf4_mul(h_block->coefficients.array[x], a);
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
        // fprintf(stderr, "flipping pos=%zu with a_max=%u and sigma_max=%ld\n", pos, a_max, sigma_max);
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
            syndrome.array[idx] ^= gf4_mul(h_block->coefficients.array[x], a_max);
            x = (0 == x) ? ctx->block_size - 1: x-1;
            ++idx;
        } while (x != h_pos);

        maybe_decoded->array[pos] ^= a_max;
    }
    gf4_array_deinit(&syndrome);
    ctx->elapsed_iterations = num_iterations;
    return false;
}
#else
bool dec_decode_symbol_flipping(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);
    long * sigmas_1 = calloc(2 * ctx->block_size, sizeof(size_t));
    long * sigmas_2 = calloc(2 * ctx->block_size, sizeof(size_t));
    long * sigmas_3 = calloc(2 * ctx->block_size, sizeof(size_t));
    long * syndrome_weights = calloc(num_iterations, sizeof(size_t));
    long * sigmas[] = {sigmas_1, sigmas_2, sigmas_3};
    char fname[100];

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long)utils_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            gf4_poly_deinit(&syndrome);
            memset(fname, 0, 100);
            sprintf(fname, "syndrome-exp_%zu.txt", ctx->index);
            FILE * outfile = fopen(fname, "w");
            for (size_t j = 0; j < i; ++j) {
                fprintf(outfile, "%zu;", syndrome_weights[j]);
            }
            fprintf(outfile, "0\n"); // correct decode means 0 syndrome weight at the end
            free(sigmas_1);
            free(sigmas_2);
            free(sigmas_3);
            free(syndrome_weights);
            fclose(outfile);
            ctx->elapsed_iterations = i;
            return true;
        }
        long sigma_max = -1;
        size_t pos = 0;
        gf4_t a_max = 0;
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
            for (gf4_t a = 1; a <= GF4_MAX_VALUE; ++a) {
                // s = s - a*h_j, where h_j is j-th column of H
                size_t x = actual_j;
                size_t idx = 0;
                long tmp_hamming_weight = 0; // hamming weight of s-a*h_j
                do {
                    gf4_t tmp = syndrome.array[idx] ^ gf4_mul(h_block->array[x], a);
                    if (0 != tmp) {
                        tmp_hamming_weight += 1;
                    }
                    x = (0 == x) ? ctx->block_size - 1: x-1;
                    ++idx;
                } while (x != actual_j);

                long sigma = syndrome_weight - tmp_hamming_weight;
                sigmas[a - 1][j] = sigma;
                if (sigma > sigma_max) {
                    sigma_max = sigma;
                    a_max = a;
                    pos = j;
                }
            }
        }
        memset(fname, 0, 100);
        sprintf(fname, "sigmas-exp_%zu-iter_%zu.txt", ctx->index, i);
        FILE * outfile = fopen(fname, "w");
        syndrome_weights[i] = (long)utils_hamming_weight(&syndrome);
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            fprintf(outfile, "%ld;", sigmas_1[j]);
        }
        fprintf(outfile, "\n");
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            fprintf(outfile, "%ld;", sigmas_2[j]);
        }
        fprintf(outfile, "\n");
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            fprintf(outfile, "%ld;", sigmas_3[j]);
        }
        fprintf(outfile, "\n");
        fclose(outfile);

        // fprintf(stderr, "flipping pos=%zu with a_max=%u and sigma_max=%ld\n", pos, a_max, sigma_max);
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
            syndrome.array[idx] ^= gf4_mul(h_block->array[x], a_max);
            x = (0 == x) ? ctx->block_size - 1: x-1;
            ++idx;
        } while (x != h_pos);

        maybe_decoded->array[pos] ^= a_max;
    }
    gf4_poly_deinit(&syndrome);

    memset(fname, 0, 100);
    sprintf(fname, "syndrome-exp_%zu.txt", ctx->index);
    FILE * outfile = fopen(fname, "w");
    for (size_t j = 0; j < num_iterations; ++j) {
        fprintf(outfile, "%zu;", syndrome_weights[j]);
    }
    fprintf(outfile, "\n");
    free(sigmas_1);
    free(sigmas_2);
    free(sigmas_3);
    free(syndrome_weights);
    fclose(outfile);

    ctx->elapsed_iterations = num_iterations;
    return false;
}
#endif