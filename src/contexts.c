#include "contexts.h"


// helper function
gf4_t poly_sum(gf4_poly_t * poly) {
    assert(NULL != poly);
    gf4_t sum = 0;
    for (size_t i = 0; i <= poly->degree; ++i) {
        sum ^= poly->coefficients[i];
    }
    return sum;
}

void contexts_init(encoding_context_t * out_enc_ctx, decoding_context_t * out_dec_ctx, size_t block_size, size_t block_weight) {
    assert(NULL != out_enc_ctx);
    assert(NULL != out_dec_ctx);
    assert(block_weight <= block_size);
    size_t capacity = block_size + 1;
    gf4_poly_t modulus = gf4_poly_init_zero(capacity);
    gf4_poly_set_coefficient(&modulus, 0, 1);
    gf4_poly_set_coefficient(&modulus, block_size, 1);

    gf4_poly_t h0 = gf4_poly_init_zero(capacity);
    gf4_poly_t h1 = gf4_poly_init_zero(capacity);
    gf4_poly_t maybe_inverse = gf4_poly_init_zero(2*capacity);
    random_weighted_gf4_poly(&h0, block_size, block_weight);
    random_weighted_gf4_poly(&h1, block_size, block_weight);

    int repetitions = 0;

    while (true) {
        while (0 == poly_sum(&h1)) {
            gf4_poly_zero_out(&h1);
            random_weighted_gf4_poly(&h1, block_size, block_weight);
        }
        bool inverted = gf4_poly_invert_slow(&maybe_inverse, &h1, &modulus);
        if (inverted) {
            gf4_poly_t tmp = gf4_poly_init_zero(2*capacity);
            gf4_poly_t div = gf4_poly_init_zero(2*capacity);
            gf4_poly_t rem = gf4_poly_init_zero(2*capacity);
            gf4_poly_mul(&tmp, &h1, &maybe_inverse);
            gf4_poly_div_rem(&div, &rem, &tmp, &modulus);
            bool correct_inverse = 0 == rem.degree && 1 == rem.coefficients[0];
            if (!correct_inverse) {
                // WTF???? this means invert function is incorrectly implemented
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
}


void contexts_load(const char * filename, size_t block_size, encoding_context_t * enc_ctx, decoding_context_t * dec_ctx) {
    assert(NULL != filename);
    assert(NULL != enc_ctx);
    assert(NULL != dec_ctx);
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
}