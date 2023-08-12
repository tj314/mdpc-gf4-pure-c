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

void dec_calculate_syndrome(gf4_vector_t *out_syndrome, gf4_vector_t *in_vector, decoding_context_t * ctx) {
    assert(NULL != out_syndrome);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(out_syndrome->capacity >= ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);
    for (size_t i = ctx->block_size; i > 0; --i) {
        gf4_t tmp = 0;
        for (size_t j = 0; j < ctx->block_size; ++j) {
            tmp ^= gf4_mul(ctx->h0.array[(i+j) % ctx->block_size], in_vector->array[j]);
            tmp ^= gf4_mul(ctx->h1.array[(i+j) % ctx->block_size], in_vector->array[ctx->block_size + j]);
        }
        out_syndrome->array[ctx->block_size - i] = tmp;
    }
    out_syndrome->length = 1;
    for (size_t i = out_syndrome->capacity; i > 0; --i) {
        if (0 != out_syndrome->array[i]) {
            out_syndrome->length = i + 1;
            break;
        }
    }
}

long dec_calculate_new_sigma(gf4_poly_t * h_block, gf4_vector_t *syndrome, long syndrome_weight, gf4_t a, size_t actual_j, decoding_context_t * ctx) {
    // s = s - a*h_j, where h_j is j-th column of H
    size_t x = actual_j;
    size_t idx = 0;
    long tmp_hamming_weight = 0; // hamming weight of s-a*h_j
    do {
        gf4_t tmp = syndrome->array[idx] ^ gf4_mul(h_block->array[x], a);
        if (0 != tmp) {
            tmp_hamming_weight += 1;
        }
        x = (0 == x) ? ctx->block_size - 1: x-1;
        ++idx;
    } while (x != actual_j);

    return syndrome_weight - tmp_hamming_weight;
}

void dec_flip_symbol(gf4_poly_t * h_block, gf4_vector_t *syndrome, gf4_vector_t *maybe_decoded, gf4_t a, size_t actual_j, size_t j, decoding_context_t * ctx) {
    size_t x = actual_j;
    size_t idx = 0;
    do {
        syndrome->array[idx] ^= gf4_mul(h_block->array[x], a);
        x = (0 == x) ? ctx->block_size - 1: x-1;
        ++idx;
    } while (x != actual_j);
    maybe_decoded->array[j] ^= a;
}

bool dec_decrypt(gf4_vector_t *out_decrypted, gf4_vector_t *in_encrypted, bool (*decode)(gf4_vector_t *, gf4_vector_t *, size_t, decoding_context_t *), size_t num_iterations, decoding_context_t * ctx) {
    bool decoding_success = decode(out_decrypted, in_encrypted, num_iterations, ctx);
    if (decoding_success) {
        memset(out_decrypted->array + ctx->block_size, 0, ctx->block_size);
        return true;
    } else {
        gf4_poly_zero_out(out_decrypted);
        return false;
    }
}