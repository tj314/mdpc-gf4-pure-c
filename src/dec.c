#include "dec.h"

void dec_calculate_syndrome(gf4_poly_t * out_syndrome, gf4_poly_t * in_vector, decoding_context_t * ctx) {
    assert(NULL != out_syndrome);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(out_syndrome->capacity >= ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);
    for (size_t i = ctx->block_size; i > 0; --i) {
        gf4_t tmp = 0;
        for (size_t j = 0; j < ctx->block_size; ++j) {
            tmp ^= gf4_mul(ctx->h0.coefficients[(i+j) % ctx->block_size], in_vector->coefficients[j]);
            tmp ^= gf4_mul(ctx->h1.coefficients[(i+j) % ctx->block_size], in_vector->coefficients[ctx->block_size + j]);
        }
        out_syndrome->coefficients[ctx->block_size - i] = tmp;
    }
    out_syndrome->degree = 0;
    for (size_t i = out_syndrome->capacity; i > 0; --i) {
        if (0 != out_syndrome->coefficients[i]) {
            out_syndrome->degree = i;
            break;
        }
    }
}


bool dec_decode_symbol_flipping(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);
#ifdef WRITE_WEIGHTS
    long * sigmas_1 = calloc(2 * ctx->block_size, sizeof(size_t));
    long * sigmas_2 = calloc(2 * ctx->block_size, sizeof(size_t));
    long * sigmas_3 = calloc(2 * ctx->block_size, sizeof(size_t));
    long * syndrome_weights = calloc(num_iterations, sizeof(size_t));
    long * sigmas[] = {sigmas_1, sigmas_2, sigmas_3};
    char fname[100];
#endif

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long)utils_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            gf4_poly_deinit(&syndrome);
#ifdef WRITE_WEIGHTS
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
#endif
            return true;
        }
        long sigma_max = -1;
        size_t pos = 0;
        gf4_t a_max = 0;
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            for (gf4_t a = 1; a <= GF4_MAX_VALUE; ++a) {
                // s = s - a*h_j, where h_j is j-th column of H
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
                long tmp_hamming_weight = 0; // hamming weight of s-a*h_j
                do {
                    gf4_t tmp = syndrome.coefficients[idx] ^ gf4_mul(h_block->coefficients[x], a);
                    if (0 != tmp) {
                        tmp_hamming_weight += 1;
                    }
                    x = (0 == x) ? ctx->block_size - 1: x-1;
                    ++idx;
                } while (x != actual_j);

                long sigma = syndrome_weight - tmp_hamming_weight;
#ifdef WRITE_WEIGHTS
                sigmas[a - 1][j] = sigma;
#endif
                if (sigma > sigma_max) {
                    sigma_max = sigma;
                    a_max = a;
                    pos = j;
                }
            }
        }
#ifdef WRITE_WEIGHTS
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
#endif
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
            syndrome.coefficients[idx] ^= gf4_mul(h_block->coefficients[x], a_max);
            x = (0 == x) ? ctx->block_size - 1: x-1;
            ++idx;
        } while (x != h_pos);

        maybe_decoded->coefficients[pos] ^= a_max;
    }
    gf4_poly_deinit(&syndrome);
#ifdef WRITE_WEIGHTS
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
#endif
    return false;
}

bool dec_decode_symbol_flipping_2(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);
#ifdef WRITE_WEIGHTS
    size_t * weights_1 = calloc(2*ctx->block_size, sizeof(size_t));
    size_t * weights_2 = calloc(2*ctx->block_size, sizeof(size_t));
    size_t * weights_3 = calloc(2*ctx->block_size, sizeof(size_t));
    size_t * syndrome_weights = calloc(num_iterations, sizeof(size_t));
    size_t * weights[] = {weights_1, weights_2, weights_3};
    char fname[100];
#endif

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long)utils_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            gf4_poly_deinit(&syndrome);
#ifdef WRITE_WEIGHTS
            memset(fname, 0, 100);
            sprintf(fname, "syndrome-exp_%zu.txt", ctx->index);
            FILE * outfile = fopen(fname, "w");
            for (size_t j = 0; j < i; ++j) {
                fprintf(outfile, "%zu;", syndrome_weights[j]);
            }
            fprintf(outfile, "0\n"); // correct decode means 0 syndrome weight at the end
            free(weights_1);
            free(weights_2);
            free(weights_3);
            free(syndrome_weights);
            fclose(outfile);
#endif
            return true;
        }
        size_t sigma_max = ctx->block_size + 1;
        size_t pos = 0;
        gf4_t a_max = 0;
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            for (gf4_t a = 1; a <= GF4_MAX_VALUE; ++a) {
                // s = s - a*h_j, where h_j is j-th column of H
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
                size_t new_syndrome_weight = 0; // hamming weight of s-a*h_j
                do {
                    gf4_t tmp = syndrome.coefficients[idx] ^ gf4_mul(h_block->coefficients[x], a);
                    if (0 != tmp) {
                        new_syndrome_weight += 1;
                    }
                    x = (0 == x) ? ctx->block_size - 1: x-1;
                    ++idx;
                } while (x != actual_j);
#ifdef WRITE_WEIGHTS
                weights[a-1][j] = new_syndrome_weight;
#endif
                if (new_syndrome_weight < sigma_max) {
                    sigma_max = new_syndrome_weight;
                    a_max = a;
                    pos = j;
                }
            }
        }
#ifdef WRITE_WEIGHTS
        memset(fname, 0, 100);
        sprintf(fname, "weights-exp_%zu-iter_%zu.txt", ctx->index, i);
        FILE * outfile = fopen(fname, "w");
        syndrome_weights[i] = utils_hamming_weight(&syndrome);
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            fprintf(outfile, "%zu;", weights_1[j]);
        }
        fprintf(outfile, "\n");
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            fprintf(outfile, "%zu;", weights_2[j]);
        }
        fprintf(outfile, "\n");
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            fprintf(outfile, "%zu;", weights_3[j]);
        }
        fprintf(outfile, "\n");
        fclose(outfile);
#endif
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
            syndrome.coefficients[idx] ^= gf4_mul(h_block->coefficients[x], a_max);
            x = (0 == x) ? ctx->block_size - 1: x-1;
            ++idx;
        } while (x != h_pos);

        maybe_decoded->coefficients[pos] ^= a_max;
    }
    gf4_poly_deinit(&syndrome);
#ifdef WRITE_WEIGHTS
    memset(fname, 0, 100);
    sprintf(fname, "syndrome-exp_%zu.txt", ctx->index);
    FILE * outfile = fopen(fname, "w");
    for (size_t j = 0; j < num_iterations; ++j) {
        fprintf(outfile, "%zu;", syndrome_weights[j]);
    }
    fprintf(outfile, "\n");
    free(weights_1);
    free(weights_2);
    free(weights_3);
    free(syndrome_weights);
    fclose(outfile);
#endif
    return false;
}

bool dec_decode_symbol_flipping_delta(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx) {
    assert(NULL != maybe_decoded);
    assert(NULL != in_vector);
    assert(NULL != ctx);
    assert(maybe_decoded->capacity >= 2*ctx->block_size);
    assert(in_vector->capacity >= 2*ctx->block_size);

    gf4_poly_t syndrome = gf4_poly_init_zero(ctx->block_size);
    dec_calculate_syndrome(&syndrome, in_vector, ctx);
    gf4_poly_copy(maybe_decoded, in_vector);

    long * sigmas = calloc(2 * ctx->block_size, sizeof(long));
    assert(NULL != sigmas);
    gf4_t * values = calloc(2 * ctx->block_size, sizeof(gf4_t));
    assert(NULL != values);

    const long DELTA = ctx->delta_setting;

    for (size_t i = 0; i < num_iterations; ++i) {
        long syndrome_weight = (long)utils_hamming_weight(&syndrome);
        if (0 == syndrome_weight) {
            free(sigmas);
            free(values);
            gf4_poly_deinit(&syndrome);
            return true;
        }
        memset(sigmas, -1, 2 * ctx->block_size);
        long sigma_max = -1;
        for (size_t j = 0; j < 2*ctx->block_size; ++j) {
            for (gf4_t a = 1; a <= GF4_MAX_VALUE; ++a) {
                // s = s - a*h_j, where h_j is j-th column of H
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
                long new_syndrome_weight = 0; // hamming weight of s-a*h_j
                do {
                    gf4_t tmp = syndrome.coefficients[idx] ^ gf4_mul(h_block->coefficients[x], a);
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
        // fprintf(stderr, "flipping pos=%zu with a_max=%u and min_syndrome_weight=%ld\n", pos, a_max, min_syndrome_weight);

        long bound = ((sigma_max - DELTA) >= 0) ? sigma_max - DELTA : 0;
        size_t flipped = 0;

        for (size_t j = 0; j < ctx->block_size; ++j) {
            if (sigmas[j] >= bound) {
                flipped += 1;
                size_t x = j;
                size_t idx = 0;
                do {
                    syndrome.coefficients[idx] ^= gf4_mul(ctx->h0.coefficients[x], values[j]);
                    x = (0 == x) ? ctx->block_size - 1: x-1;
                    ++idx;
                } while (x != j);
                maybe_decoded->coefficients[j] ^= values[j];
            }
        }
        for (size_t j = 0; j < ctx->block_size; ++j) {
            if (sigmas[ctx->block_size + j] >= bound) {
                flipped += 1;
                size_t x = j;
                size_t idx = 0;
                do {
                    syndrome.coefficients[idx] ^= gf4_mul(ctx->h1.coefficients[x], values[ctx->block_size + j]);
                    x = (0 == x) ? ctx->block_size - 1: x-1;
                    ++idx;
                } while (x != j);
                maybe_decoded->coefficients[ctx->block_size + j] ^= values[ctx->block_size + j];
            }
        }
        fprintf(stderr, "flipped: %zu\n", flipped);
    }
    free(sigmas);
    free(values);
    gf4_poly_deinit(&syndrome);
    return false;
}

bool dec_decrypt(gf4_poly_t * out_decrypted, gf4_poly_t * in_encrypted, bool (*decode)(gf4_poly_t*, gf4_poly_t *, size_t, decoding_context_t *), size_t num_iterations, decoding_context_t * ctx) {
    bool decoding_success = decode(out_decrypted, in_encrypted, num_iterations, ctx);
    if (decoding_success) {
        memset(out_decrypted->coefficients + ctx->block_size, 0, ctx->block_size);
        return true;
    } else {
        gf4_poly_zero_out(out_decrypted);
        return false;
    }
}