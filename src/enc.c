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


void enc_encode(gf4_poly_t * out_encoded, gf4_poly_t * in_message, encoding_context_t * ctx) {
    assert(NULL != out_encoded);
    assert(NULL != in_message);
    assert(NULL != ctx);
    assert(out_encoded->capacity >= 2*ctx->block_size);
    assert(in_message->capacity >= ctx->block_size);
    memcpy(out_encoded->coefficients, in_message->coefficients, ctx->block_size);
    for (size_t i = ctx->block_size; i > 0; --i) {
        gf4_t tmp = 0;
        for (size_t j = 0; j < ctx->block_size; ++j) {
            tmp ^= gf4_mul(in_message->coefficients[j], ctx->second_block_G.coefficients[(i + j) % ctx->block_size]);
        }
        out_encoded->coefficients[2*ctx->block_size - i] = tmp;
    }
    gf4_poly_adjust_degree(out_encoded, out_encoded->capacity - 1);
}

void enc_encrypt(gf4_poly_t * out_encrypted, gf4_poly_t * in_message, size_t num_errors, encoding_context_t * ctx) {
    enc_encode(out_encrypted, in_message, ctx);
    gf4_poly_t err = gf4_poly_init_zero(2*ctx->block_size);
    random_weighted_gf4_poly(&err, 2*ctx->block_size, num_errors);
#ifdef WRITE_WEIGHTS
    char fname[100] = {0};
    sprintf(fname, "errorvec-exp_%zu.txt", ctx->index);
    FILE * outfile = fopen(fname, "w");
    for (size_t i = 0; i < 2*ctx->block_size; ++i) {
        fprintf(outfile, "%u;", err.coefficients[i]);
    }
    fprintf(outfile, "\n");
    fclose(outfile);
#endif
    gf4_poly_add_inplace(out_encrypted, &err);
    gf4_poly_deinit(&err);
}