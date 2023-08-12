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

#include "contexts.h"
#include "utils.h"


// helper function
gf4_t contexts_poly_sum(gf4_poly_t * poly) {
    assert(NULL != poly);
    gf4_t sum = 0;
    for (size_t i = 0; i < poly->length; ++i) {
        sum ^= poly->array[i];
    }
    return sum;
}

void contexts_init(encoding_context_t * out_enc_ctx, decoding_context_t * out_dec_ctx, size_t block_size, size_t block_weight) {
    assert(NULL != out_enc_ctx);
    assert(NULL != out_dec_ctx);
    assert(block_weight <= block_size);

    // settings not required by all decoders
    out_dec_ctx->threshold = NULL;
    out_dec_ctx->elapsed_iterations = 0;
    out_dec_ctx->delta_setting = -1;

    // generate keys
    size_t capacity = block_size + 1;
    gf4_poly_t modulus = gf4_poly_init_zero(capacity);
    gf4_poly_set_coefficient(&modulus, 0, 1);
    gf4_poly_set_coefficient(&modulus, block_size, 1);

    gf4_poly_t h0 = gf4_poly_init_zero(capacity);
    gf4_poly_t h1 = gf4_poly_init_zero(capacity);
    gf4_poly_t maybe_inverse = gf4_poly_init_zero(2*capacity);
    random_weighted_gf4_vector(&h0, block_size, block_weight);
    random_weighted_gf4_vector(&h1, block_size, block_weight);

    int repetitions = 0;

    while (true) {
        while (0 == contexts_poly_sum(&h1)) {
            gf4_poly_zero_out(&h1);
            random_weighted_gf4_vector(&h1, block_size, block_weight);
        }
        bool inverted = gf4_poly_invert_slow(&maybe_inverse, &h1, &modulus);
        if (inverted) {
            gf4_poly_t tmp = gf4_poly_init_zero(2*capacity);
            gf4_poly_t div = gf4_poly_init_zero(2*capacity);
            gf4_poly_t rem = gf4_poly_init_zero(2*capacity);
            gf4_poly_mul(&tmp, &h1, &maybe_inverse);
            gf4_poly_div_rem(&div, &rem, &tmp, &modulus);
            bool correct_inverse = 0 == gf4_poly_get_degree(&rem) && 1 == rem.array[0];
            if (!correct_inverse) {
                // WTF???? this means invert function is incorrectly implemented
                fprintf(stderr, "contexts_init: WTF? invert function is incorrectly implemented!\n");
                gf4_poly_deinit(&tmp);
                gf4_poly_deinit(&rem);
                gf4_poly_deinit(&div);
                gf4_poly_deinit(&modulus);
                gf4_poly_deinit(&maybe_inverse);
                gf4_poly_deinit(&h0);
                gf4_poly_deinit(&h1);
                exit(-1);
            }

            // second_block_G_poly = (h0_poly * inverse) % modulus ;
            gf4_poly_zero_out(&tmp);
            gf4_poly_zero_out(&rem);
            gf4_poly_zero_out(&div);
            gf4_poly_mul(&tmp, &h0, &maybe_inverse);
            gf4_poly_div_rem(&div, &rem, &tmp, &modulus);
            gf4_poly_deinit(&tmp);
            gf4_poly_deinit(&div);
            out_enc_ctx->block_size = block_size;
            out_enc_ctx->second_block_G = rem;
            out_dec_ctx->block_size = block_size;
            out_dec_ctx->h0 = h0;
            out_dec_ctx->h1 = h1;
            gf4_poly_deinit(&modulus);
            gf4_poly_deinit(&maybe_inverse);
            return;
        }
        gf4_poly_zero_out(&maybe_inverse);
        gf4_poly_zero_out(&h1);
        repetitions++;
        if (repetitions > 20) {
            random_force_reseed();
            repetitions = 0;
        }
    }
}

void contexts_deinit(encoding_context_t * enc_ctx, decoding_context_t * dec_ctx) {
    assert(NULL != enc_ctx);
    assert(NULL != dec_ctx);
    gf4_poly_deinit(&(enc_ctx->second_block_G));
    gf4_poly_deinit(&(dec_ctx->h0));
    gf4_poly_deinit(&(dec_ctx->h1));
    enc_ctx->block_size = 0;
    dec_ctx->block_size = 0;
}

void contexts_save(const char * filename, encoding_context_t * enc_ctx, decoding_context_t * dec_ctx) {
    assert(NULL != filename);
    assert(NULL != enc_ctx);
    assert(NULL != dec_ctx);
    /*
    FILE * output = fopen(filename, "wb+");
    if (NULL == output) {
        fprintf(stderr, "contexts_save: Output file couldn't be created!\n");
        exit(-1);
    }
    size_t num_nonzero = utils_hamming_weight(&enc_ctx->second_block_G);
    fwrite(&num_nonzero, sizeof(size_t), 1, output);
    for (size_t i = 0; i <= enc_ctx->second_block_G.degree; ++i) {
        if (0 != enc_ctx->second_block_G.coefficients[i]) {
            fwrite(&i, sizeof(size_t), 1, output);
            fwrite(&enc_ctx->second_block_G.coefficients[i], sizeof(uint8_t), 1, output);
        }
    }

    num_nonzero = utils_hamming_weight(&dec_ctx->h0);
    fwrite(&num_nonzero, sizeof(size_t), 1, output);
    for (size_t i = 0; i <= dec_ctx->h0.degree; ++i) {
        if (0 != dec_ctx->h0.coefficients[i]) {
            fwrite(&i, sizeof(size_t), 1, output);
            fwrite(&dec_ctx->h0.coefficients[i], sizeof(uint8_t), 1, output);
        }
    }

    num_nonzero = utils_hamming_weight(&dec_ctx->h1);
    fwrite(&num_nonzero, sizeof(size_t), 1, output);
    for (size_t i = 0; i <= dec_ctx->h1.degree; ++i) {
        if (0 != dec_ctx->h1.coefficients[i]) {
            fwrite(&i, sizeof(size_t), 1, output);
            fwrite(&dec_ctx->h1.coefficients[i], sizeof(uint8_t), 1, output);
        }
    }

    fclose(output);
    */
    FILE * output = fopen(filename, "w+");
    if (NULL == output) {
        fprintf(stderr, "contexts_save: Output file couldn't be created!\n");
        exit(-1);
    }
    size_t num_nonzero = utils_hamming_weight(&enc_ctx->second_block_G);
    fprintf(output, "%zu\n", num_nonzero);
    for (size_t i = 0; i < enc_ctx->second_block_G.length; ++i) {
        if (0 != enc_ctx->second_block_G.array[i]) {
            fprintf(output, "%zu %hhu\n", i, enc_ctx->second_block_G.array[i]);
        }
    }

    num_nonzero = utils_hamming_weight(&dec_ctx->h0);
    fprintf(output, "%zu\n", num_nonzero);
    for (size_t i = 0; i < dec_ctx->h0.length; ++i) {
        if (0 != dec_ctx->h0.array[i]) {
            fprintf(output, "%zu %hhu\n", i, dec_ctx->h0.array[i]);
        }
    }

    num_nonzero = utils_hamming_weight(&dec_ctx->h1);
    fprintf(output, "%zu\n", num_nonzero);
    for (size_t i = 0; i < dec_ctx->h1.length; ++i) {
        if (0 != dec_ctx->h1.array[i]) {
            fprintf(output, "%zu %hhu\n", i, dec_ctx->h1.array[i]);
        }
    }

    fclose(output);
}


void contexts_load(const char * filename, size_t block_size, encoding_context_t * enc_ctx, decoding_context_t * dec_ctx) {
    assert(NULL != filename);
    assert(NULL != enc_ctx);
    assert(NULL != dec_ctx);
    /*
    FILE * input = fopen(filename, "rb");
    if (NULL == input) {
        fprintf(stderr, "contexts_load: Input file doesn't exist!\n");
        exit(-1);
    }

    enc_ctx->block_size = block_size;
    dec_ctx->block_size = block_size;
    enc_ctx->second_block_G = gf4_poly_init_zero(block_size);
    dec_ctx->h0 = gf4_poly_init_zero(block_size);
    dec_ctx->h1 = gf4_poly_init_zero(block_size);

    size_t deg, num_nonzero;
    uint8_t coeff;

    fread(&num_nonzero, sizeof(size_t), 1, input);
    for (size_t i = 0; i < num_nonzero; ++i) {
        fread(&deg, sizeof(size_t), 1, input);
        fread(&coeff, sizeof(uint8_t), 1, input);
        enc_ctx->second_block_G.coefficients[deg] = coeff;
        enc_ctx->second_block_G.degree = deg;
    }

    fread(&num_nonzero, sizeof(size_t), 1, input);
    for (size_t i = 0; i < num_nonzero; ++i) {
        fread(&deg, sizeof(size_t), 1, input);
        fread(&coeff, sizeof(uint8_t), 1, input);
        dec_ctx->h0.coefficients[deg] = coeff;
        dec_ctx->h0.degree = deg;
    }

    fread(&num_nonzero, sizeof(size_t), 1, input);
    for (size_t i = 0; i < num_nonzero; ++i) {
        fread(&deg, sizeof(size_t), 1, input);
        fread(&coeff, sizeof(uint8_t), 1, input);
        dec_ctx->h1.coefficients[deg] = coeff;
        dec_ctx->h1.degree = deg;
    }

    fclose(input);
    */
    FILE * input = fopen(filename, "r");
    if (NULL == input) {
        fprintf(stderr, "contexts_load: Input file doesn't exist!\n");
        exit(-1);
    }

    enc_ctx->block_size = block_size;
    dec_ctx->block_size = block_size;
    enc_ctx->second_block_G = gf4_poly_init_zero(block_size);
    dec_ctx->h0 = gf4_poly_init_zero(block_size);
    dec_ctx->h1 = gf4_poly_init_zero(block_size);

    size_t deg, num_nonzero;
    gf4_t coeff;

    fscanf(input, "%zu", &num_nonzero);
    fgetc(input); // remove line ending
    for (size_t i = 0; i < num_nonzero; ++i) {
        fscanf(input, "%zu", &deg);
        fscanf(input, "%hhu", &coeff);
        fgetc(input); // remove line ending
        enc_ctx->second_block_G.array[deg] = coeff;
        enc_ctx->second_block_G.length = deg + 1;
    }

    fscanf(input, "%zu", &num_nonzero);
    fgetc(input); // remove line ending
    for (size_t i = 0; i < num_nonzero; ++i) {
        fscanf(input, "%zu", &deg);
        fscanf(input, "%hhu", &coeff);
        fgetc(input); // remove line ending
        dec_ctx->h0.array[deg] = coeff;
        dec_ctx->h0.length = deg + 1;
    }

    fscanf(input, "%zu", &num_nonzero);
    fgetc(input); // remove line ending
    for (size_t i = 0; i < num_nonzero; ++i) {
        fscanf(input, "%zu", &deg);
        fscanf(input, "%hhu", &coeff);
        fgetc(input); // remove line ending
        dec_ctx->h1.array[deg] = coeff;
        dec_ctx->h1.length = deg + 1;
    }

    fclose(input);
}