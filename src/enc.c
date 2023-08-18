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

#include "enc.h"


void enc_encode(gf4_array_t *out_encoded, gf4_array_t *in_message, encoding_context_t * ctx) {
    assert(NULL != out_encoded);
    assert(NULL != in_message);
    assert(NULL != ctx);
    assert(out_encoded->capacity >= 2*ctx->block_size);
    assert(in_message->capacity >= ctx->block_size);
    memcpy(out_encoded->array, in_message->array, ctx->block_size);
    for (size_t i = ctx->block_size; i > 0; --i) {
        gf4_t tmp = 0;
        for (size_t j = 0; j < ctx->block_size; ++j) {
            tmp ^= gf4_mul(in_message->array[j], ctx->second_block_G.coefficients.array[(i + j) % ctx->block_size]);
        }
        out_encoded->array[2*ctx->block_size - i] = tmp;
    }
}

void enc_encrypt(gf4_array_t *out_encrypted, gf4_array_t *in_message, size_t num_errors, encoding_context_t * ctx) {
    assert(NULL != out_encrypted);
    assert(NULL != in_message);
    assert(NULL != ctx);
    assert(out_encrypted->capacity >= 2*ctx->block_size);
    assert(in_message->capacity >= ctx->block_size);
    enc_encode(out_encrypted, in_message, ctx);
    gf4_array_t err = gf4_array_init(2*ctx->block_size, true);
    random_weighted_gf4_array(&err, 2 * ctx->block_size, num_errors);
#ifdef WRITE_WEIGHTS
    char fname[100] = {0};
    sprintf(fname, "errorvec-exp_%zu.txt", ctx->index);
    FILE * outfile = fopen(fname, "w");
    for (size_t i = 0; i < 2*ctx->block_size; ++i) {
        fprintf(outfile, "%u;", err.coefficients.array[i]);
    }
    fprintf(outfile, "\n");
    fclose(outfile);
#endif
    for (size_t i = 0; i < 2*ctx->block_size; ++i) {
        out_encrypted->array[i] = gf4_add(out_encrypted->array[i], err.array[i]);
    }
    gf4_array_deinit(&err);
}