#ifndef RUNTESTS
#include <stdio.h>
#include <pthread.h>
#include "src/gf4.h"
#include "src/gf4_poly.h"
#include "src/contexts.h"
#include "src/enc.h"
#include "src/dec.h"

void experiment0(size_t block_size, size_t num_errors, size_t M, const char * out_filename) {
    encoding_context_t ec;
    decoding_context_t dc;
    contexts_load("keys.bin", block_size, &ec, &dc);

    size_t bound = (size_t)(block_size / 2);

    gf4_poly_t syndrome = gf4_poly_init_zero(block_size);
    gf4_poly_t err_vector = gf4_poly_init_zero(2*block_size);

    FILE * output = fopen(out_filename, "w+");

    for (size_t d = 1; d <= bound; ++d) {
        fprintf(output, "\ndistance: %zu\n", d);
        fprintf(output, "===================\n\n");
        for (size_t i = 0; i < M; ++i) {
            random_weighted_gf4_poly_pairs_of_ones(&err_vector, block_size, num_errors, d);
            dec_calculate_syndrome(&syndrome, &err_vector, &dc);
            size_t syndrome_weight = utils_hamming_weight(&syndrome);

            fprintf(output, "err: ");
            gf4_poly_coeff_print(&err_vector, 2*block_size, output, "\n");
            fprintf(output, "syn: ");
            gf4_poly_coeff_print(&syndrome, block_size, output, "\n");
            fprintf(output, "syndrome weight: %zu\n\n", syndrome_weight);

            gf4_poly_zero_out(&err_vector);
            gf4_poly_zero_out(&syndrome);
        }
    }

    fclose(output);

    gf4_poly_deinit(&syndrome);
    gf4_poly_deinit(&err_vector);
    contexts_deinit(&ec, &dc);
}

void experiment1(size_t block_size, size_t num_errors, size_t M, const char * out_filename) {
    encoding_context_t ec;
    decoding_context_t dc;
    contexts_load("keys.bin", block_size, &ec, &dc);

    size_t bound = (size_t)(block_size / 2);

    gf4_poly_t syndrome = gf4_poly_init_zero(block_size);
    gf4_poly_t err_vector = gf4_poly_init_zero(2*block_size);

    FILE * output = fopen(out_filename, "w+");

    for (size_t d = 1; d <= bound; ++d) {
        fprintf(output, "\ndistance: %zu\n", d);
        fprintf(output, "===================\n\n");
        for (size_t i = 0; i < M; ++i) {
            random_weighted_gf4_poly_pairs_of_one_alpha(&err_vector, block_size, num_errors, d);
            dec_calculate_syndrome(&syndrome, &err_vector, &dc);
            size_t syndrome_weight = utils_hamming_weight(&syndrome);

            fprintf(output, "err: ");
            gf4_poly_coeff_print(&err_vector, 2*block_size, output, "\n");
            fprintf(output, "syn: ");
            gf4_poly_coeff_print(&syndrome, block_size, output, "\n");
            fprintf(output, "syndrome weight: %zu\n\n", syndrome_weight);

            gf4_poly_zero_out(&err_vector);
            gf4_poly_zero_out(&syndrome);
        }
    }

    fclose(output);

    gf4_poly_deinit(&syndrome);
    gf4_poly_deinit(&err_vector);
    contexts_deinit(&ec, &dc);
}

void gen_keys(size_t block_size, size_t block_weight) {
    encoding_context_t ec;
    decoding_context_t dc;
    contexts_init(&ec, &dc, block_size, block_weight);
    contexts_save("keys.bin", &ec, &dc);
    contexts_deinit(&ec, &dc);
}


void run_tests(FILE * file, size_t num_keys, size_t num_messages, size_t block_size, size_t block_weight, size_t num_errors, size_t num_iterations, size_t decoder, size_t delta) {
    /*
    size_t block_size = 2339;
    size_t block_weight = 37;
    size_t num_errors = 84;
    size_t num_iterations = 150;
    */
    bool (*decode_function)(gf4_poly_t*, gf4_poly_t *, size_t, decoding_context_t *);
    switch (decoder) {
        case 0:
            decode_function = &dec_decode_symbol_flipping;
            break;
        case 1:
            decode_function = &dec_decode_symbol_flipping_2;
            break;
        case 2:
            decode_function = &dec_decode_symbol_flipping_delta;
            break;
        default:
            decode_function = NULL;
            exit(-1);
    }
    gf4_poly_t plaintext = gf4_poly_init_zero(block_size);
    gf4_poly_t ciphertext = gf4_poly_init_zero(2 * block_size);
    gf4_poly_t decrypted = gf4_poly_init_zero(2 * block_size);
    size_t failure_counter = 0;
    for (size_t i = 0; i < num_keys; ++i) {
        fprintf(stderr, "regen keys!\n");
        encoding_context_t ec;
        decoding_context_t dc;
        contexts_init(&ec, &dc, block_size, block_weight);

        dc.delta_setting = (long)delta; // this setting does not matter for decoders that do not use delta // TODO: remove

        for (size_t j = 0; j < num_messages; ++j) {
            random_gf4_poly(&plaintext, block_size);
            enc_encrypt(&ciphertext, &plaintext, num_errors, &ec);
            bool decryption_success = dec_decrypt(&decrypted, &ciphertext, dec_decode_symbol_flipping, num_iterations,
                                                  &dc);
            if (decryption_success) {
                fprintf(stderr, "progress: %zu / %zu,  SUCCESS\n", j + 1, num_messages*num_keys);
            } else {
                failure_counter += 1;
                fprintf(stderr, "progress: %zu / %zu,  FAILURE\n", j + 1, num_messages*num_keys);
            }
            gf4_poly_zero_out(&plaintext);
            gf4_poly_zero_out(&ciphertext);
            gf4_poly_zero_out(&decrypted);
        }
        contexts_deinit(&ec, &dc);
    }
    gf4_poly_deinit(&plaintext);
    gf4_poly_deinit(&ciphertext);
    gf4_poly_deinit(&decrypted);
    fprintf(file, "num failures: %zu / %zu\n", failure_counter, num_messages);
}


void print_usage() {
    fprintf(stderr, "Usage: ./mdpc-gf4 NUM_KEYS NUM_MSGS BLOCK_SIZE BLOCK_WEIGHT NUM_ERRORS NUM_ITERS DECODER [DELTA]\n");
    fprintf(stderr, "e.g.:  ./mdpc-gf4 10 100 2293 37 88 200 0\n");
    fprintf(stderr, "\nParameter values:\n");
    fprintf(stderr, "NUM_KEYS:     positive integer, number of key pairs to generate\n");
    fprintf(stderr, "NUM_MSGS:     positive integer, number of messages to generate, encrypt and decrypt for every key pair\n");
    fprintf(stderr, "BLOCK_SIZE:   positive integer, size of cyclic matrix block, 2293 or 2339 is recommended\n");
    fprintf(stderr, "BLOCK_WEIGHT: positive integer, hamming weight of cyclic matrix block, 37 is recommended\n");
    fprintf(stderr, "NUM_ERRORS:   positive integer, number of errors in the error vector\n");
    fprintf(stderr, "NUM_ITERS:    positive integer, number of decoding iterations\n");
    fprintf(stderr, "DECODER:      positive integer, one of the options listed bellow\n");
    fprintf(stderr, "DELTA:        optional positive integer, used as delta with decoder 2\n");
    fprintf(stderr, "\nPossible decoders:\n");
    fprintf(stderr, "0 --> SF v1\n");
    fprintf(stderr, "1 --> SF v2\n");
    fprintf(stderr, "2 --> SF with delta\n");
    fprintf(stderr, "3 --> SF with thr\n");
}

int main(int argc, char ** argv) {
    /*
    size_t block_size = 2293;
    size_t block_weight = 37;
    size_t num_iterations = 150;
    size_t num_messages = 200;
    */
    if (8 != argc && 9 != argc) {
        print_usage();
        return 0;
    }
    size_t num_keys;
    size_t num_messages;
    size_t block_size;
    size_t block_weight;
    size_t num_errors;
    size_t num_iterations;
    size_t decoder;
    size_t delta;
    num_keys = atol(argv[1]);
    num_messages = atol(argv[2]);
    block_size = atol(argv[3]);
    block_weight = atol(argv[4]);
    num_errors = atol(argv[5]);
    num_iterations = atol(argv[6]);
    decoder = atol(argv[7]);
    fprintf(stderr, "Params are: %zu %zu %zu %zu %zu %zu %zu", num_keys, num_messages, block_size, block_weight, num_errors, num_iterations, decoder);
    if (8 == argc) {
        fprintf(stderr, "\n");
        run_tests(stderr, num_keys, num_messages, block_size, block_weight, num_errors, num_iterations, decoder, 0);
    } else { // delta is provided
        delta = atol(argv[8]);
        fprintf(stderr, " %zu\n", delta);
        run_tests(stderr, num_keys, num_messages, block_size, block_weight, num_errors, num_iterations, decoder, delta);
    }
    return 0;
}
#else
#include <stdio.h>
#include "src/tests.h"

int main() {
    test_dec_calculate_syndrome();
    return 0;
}
#endif