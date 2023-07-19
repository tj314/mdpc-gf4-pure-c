#define GJS

#ifdef RUNTESTS // run unit tests defined in test.h (new tests must be called here manually!)
#include <stdio.h>
#include "src/tests.h"

int main() {
    test_gf4_is_in_range();
    test_gf4_to_str();
    test_gf4_add();
    test_gf4_mul();
    test_gf4_div();
    test_utils_hamming_weight();
    test_gf4_poly_cyclic_shift_right_inplace();
    test_gf4_square_matrix_make_cyclic_matrix();
    test_contexts_init();
    test_contexts_save_load();
    test_enc_encode();
    test_enc_encrypt();
    test_dec_calculate_syndrome();
    fprintf(stderr, "All tests passed!\n");
    return 0;
}

#elif defined(WRITE_WIGHTS) // write syndrome weights and sigmas to files for later analysis
void test_syndromes() {
    // define WRITE_WEIGHTS
    size_t block_size = 2339;
    size_t block_weight = 37;
    size_t num_iterations = 150;
    size_t num_errors = 84;
    gf4_poly_t plaintext = gf4_poly_init_zero(block_size);
    gf4_poly_t ciphertext = gf4_poly_init_zero(2 * block_size);
    gf4_poly_t decrypted = gf4_poly_init_zero(2 * block_size);
    size_t counter = 0;
    for (size_t key = 0; key < 10; ++key) {
        encoding_context_t ec;
        decoding_context_t dc;
        contexts_init(&ec, &dc, block_size, block_weight);
        for (size_t msg = 0; msg < 10; ++msg) {
            fprintf(stderr, "Progress: %zu/100\n", counter);
            dc.index = counter;
            ec.index = counter;
            random_gf4_poly(&plaintext, block_size);
            enc_encrypt(&ciphertext, &plaintext, num_errors, &ec);
            bool success = dec_decrypt(&decrypted, &ciphertext, dec_decode_symbol_flipping, num_iterations,&dc);
            if (success) {
                counter++;
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
}

int main() {
    test_syndromes();
    return 0;
}

#elif defined(TEST_ITERATIONS) // find number of decoder iterations

#include "src/dec.h"
#include "src/enc.h"

void test_iterations(size_t decoder, size_t opt) {
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
        case 3:
            decode_function = &dec_decode_symbol_flipping_threshold;
            break;
        case 4:
            decode_function = &dec_decode_symbol_flipping_bg;
            break;
        default:
            decode_function = NULL;
            exit(-1);
    }

    // select threshold if relevant
    long (*threshold_function)(long) = NULL;
    if (3 == decoder) {
        switch (opt) {
            case 0:
                threshold_function = &dec_calculate_threshold_0;
                break;
            case 1:
                threshold_function = &dec_calculate_threshold_1;
                break;
            case 2:
                threshold_function = &dec_calculate_threshold_2;
                break;
            case 3:
                threshold_function = &dec_calculate_threshold_3;
                break;
            case 4:
                threshold_function = &dec_calculate_threshold_4;
                break;
            case 5:
                threshold_function = &dec_calculate_threshold_5;
                break;
            default:
                fprintf(stderr, "ERROR: incorrect OPT value!\n");
                exit(-1);
        }
    }

    if (4 == decoder) {
        threshold_function = &dec_calculate_threshold_3;
    }

    size_t num_keys = 20;
    size_t num_messages = 100;
    size_t block_size = 2339;
    size_t block_weight = 37;
    size_t num_errors = 84;
    size_t num_iterations = 200;

    gf4_poly_t plaintext = gf4_poly_init_zero(block_size);
    gf4_poly_t ciphertext = gf4_poly_init_zero(2 * block_size);
    gf4_poly_t decrypted = gf4_poly_init_zero(2 * block_size);
    size_t * iteration_counter = malloc(num_keys*num_messages* sizeof(size_t));
    size_t pos = 0;
    for (size_t i = 0; i < num_keys; ++i) {
        fprintf(stderr, "regen keys!\n");
        encoding_context_t ec;
        decoding_context_t dc;
        contexts_init(&ec, &dc, block_size, block_weight);

        dc.delta_setting = (long)opt;
        dc.threshold = threshold_function;

        for (size_t j = 0; j < num_messages; ++j) {
            RERUN:
            random_gf4_poly(&plaintext, block_size);
            enc_encrypt(&ciphertext, &plaintext, num_errors, &ec);
            bool decryption_success = dec_decrypt(&decrypted, &ciphertext, decode_function, num_iterations,
                                                  &dc);
            if (decryption_success) {
                fprintf(stderr, "progress: %zu / %zu,  SUCCESS\n", pos+1, num_messages*num_keys);
                gf4_poly_zero_out(&plaintext);
                gf4_poly_zero_out(&ciphertext);
                gf4_poly_zero_out(&decrypted);
                iteration_counter[pos] = dc.elapsed_iterations;
                pos++;
            } else {
                gf4_poly_zero_out(&plaintext);
                gf4_poly_zero_out(&ciphertext);
                gf4_poly_zero_out(&decrypted);
                fprintf(stderr, "progress: %zu / %zu,  FAILURE, RERUN\n", pos+1, num_messages*num_keys);
                goto RERUN;
            }
        }
        contexts_deinit(&ec, &dc);
    }
    gf4_poly_deinit(&plaintext);
    gf4_poly_deinit(&ciphertext);
    gf4_poly_deinit(&decrypted);

    char fname[100] = {0};
    if (2 == decoder) {
        sprintf(fname, "iteracie-dec_%zu-delta_%zu.txt", decoder, opt);
    } else if (3 == decoder) {
        sprintf(fname, "iteracie-dec_%zu-threshold_%zu.txt", decoder, opt);
    } else {
        sprintf(fname, "iteracie-dec_%zu-opt_%zu.txt", decoder, opt);
    }

    FILE * out = fopen(fname, "w");
    for (size_t i = 0; i < pos; ++i) {
        fprintf(out, "%zu;", iteration_counter[i]);
    }
    fprintf(out, "\n");
    fclose(out);
    free(iteration_counter);
}

int main(int argc, char ** argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: ./mdpc-gf4 DECODER OPT\n");
        fprintf(stderr, "DECODER: 0 --> SF\n");
        fprintf(stderr, "         2 --> SF with DELTA\n");
        fprintf(stderr, "         3 --> SF with threshold\n");
        fprintf(stderr, "OPT setting is required, but is ignored in decoders that do not need it.\n");
        fprintf(stderr, "With decoder 2, OPT is used as delta setting\n");
        fprintf(stderr, "With decoder 3, OPT is used as index of threshold function\n");
        return -1;
    }
    size_t decoder = atol(argv[1]);
    size_t delta = atol(argv[2]);
    fprintf(stderr, "Settings are: %zu %zu\n", decoder, delta);
    test_iterations(decoder, delta);
}

#elif defined(GJS)  // TODO, for the love of god do not run this, bad things will happen, monsters will crawl from under your bed

#include <stdio.h>
#include "src/gf4_poly.h"
#include "src/contexts.h"
#include "src/dec.h"


int main() {
    size_t ** mults_same = malloc(4* sizeof(size_t *));
    assert(NULL != mults_same);
    size_t ** mults_diff = malloc(4* sizeof(size_t *));
    assert(NULL != mults_diff);
    for (size_t i = 0; i < 4; ++i) {
        mults_same[i] = calloc(2339, sizeof(size_t));
        assert(NULL != mults_same[i]);
        mults_diff[i] = calloc(2339, sizeof(size_t));
        assert(NULL != mults_diff[i]);
    }
    encoding_context_t ec;
    decoding_context_t dc;
    contexts_init(&ec, &dc, 2339, 37);

    utils_get_distance_multiplicities_h0(mults_same, mults_diff, &dc);

    fprintf(stderr, "Mults same symbols:\n");
    for (size_t i = 0; i < 4; ++i) {
        fprintf(stderr, "mult = %zu: ", i);
        for (size_t j = 0; j < 20; ++j) {
            fprintf(stderr, "%zu, ", mults_same[i][j]);
        }
        fprintf(stderr, "%zu\n", mults_same[i][20]);
    }

    fprintf(stderr, "Mults diff symbols:\n");
    for (size_t i = 0; i < 4; ++i) {
        fprintf(stderr, "mult = %zu: ", i);
        for (size_t j = 0; j < 20; ++j) {
            fprintf(stderr, "%zu, ", mults_diff[i][j]);
        }
        fprintf(stderr, "%zu\n", mults_diff[i][20]);
    }

    for (size_t i = 0; i < 4; ++i) {
        free(mults_same[i]);
        free(mults_diff[i]);
    }
    free(mults_same);
    free(mults_diff);
    contexts_deinit(&ec, &dc);
    return 0;
}

#else // DFR tests

#include <stdio.h>
#include "src/gf4.h"
#include "src/gf4_poly.h"
#include "src/contexts.h"
#include "src/enc.h"
#include "src/dec.h"

void run_tests(FILE * file, size_t num_keys, size_t num_messages, size_t block_size, size_t block_weight, size_t num_errors, size_t num_iterations, size_t decoder, size_t opt) {
    // select decoder
    bool (*decode_function)(gf4_poly_t*, gf4_poly_t *, size_t, decoding_context_t *) = NULL;
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
        case 3:
            decode_function = &dec_decode_symbol_flipping_threshold;
            break;
        case 4:
            decode_function = &dec_decode_symbol_flipping_bg;
            break;
        default:
            decode_function = NULL;
            exit(-1);
    }

    // select threshold if relevant
    // select threshold if relevant
    long (*threshold_function)(long) = NULL;
    if (3 == decoder) {
        switch (opt) {
            case 0:
                threshold_function = &dec_calculate_threshold_0;
                break;
            case 1:
                threshold_function = &dec_calculate_threshold_1;
                break;
            case 2:
                threshold_function = &dec_calculate_threshold_2;
                break;
            case 3:
                threshold_function = &dec_calculate_threshold_3;
                break;
            case 4:
                threshold_function = &dec_calculate_threshold_4;
                break;
            case 5:
                threshold_function = &dec_calculate_threshold_5;
                break;
            default:
                fprintf(stderr, "ERROR: incorrect OPT value!\n");
                exit(-1);
        }
    }

    if (4 == decoder) {
        threshold_function = &dec_calculate_threshold_3;
    }

    // tests
    gf4_poly_t plaintext = gf4_poly_init_zero(block_size);
    gf4_poly_t ciphertext = gf4_poly_init_zero(2 * block_size);
    gf4_poly_t decrypted = gf4_poly_init_zero(2 * block_size);
    size_t failure_counter = 0;
    size_t pos = 0;
    for (size_t i = 0; i < num_keys; ++i) {
        fprintf(stderr, "regen keys! current num failures: %zu / %zu\n", failure_counter, num_keys*num_messages);
        encoding_context_t ec;
        decoding_context_t dc;
        contexts_init(&ec, &dc, block_size, block_weight);

        if (3 == decoder || 4 == decoder) {
            dc.threshold = threshold_function;
            if (4 == decoder) {
                dc.delta_setting = (long)opt;
            }
        } else {
            dc.delta_setting = (long)opt;
        }

        for (size_t j = 0; j < num_messages; ++j) {
            random_gf4_poly(&plaintext, block_size);
            enc_encrypt(&ciphertext, &plaintext, num_errors, &ec);
            bool decryption_success = dec_decrypt(&decrypted, &ciphertext, decode_function, num_iterations,
                                                  &dc);
            if (decryption_success) {
                fprintf(stderr, "progress: %zu / %zu,  SUCCESS\n", pos+1, num_messages*num_keys);
            } else {
                failure_counter += 1;
                fprintf(stderr, "progress: %zu / %zu,  FAILURE\n", pos+1, num_messages*num_keys);
            }
            gf4_poly_zero_out(&plaintext);
            gf4_poly_zero_out(&ciphertext);
            gf4_poly_zero_out(&decrypted);
            pos++;
        }
        contexts_deinit(&ec, &dc);
    }
    gf4_poly_deinit(&plaintext);
    gf4_poly_deinit(&ciphertext);
    gf4_poly_deinit(&decrypted);
    fprintf(file, "num failures: %zu / %zu\n", failure_counter, num_keys*num_messages);
}

void print_usage() {
    fprintf(stderr, "Usage: ./mdpc-gf4 NUM_KEYS NUM_MSGS BLOCK_SIZE BLOCK_WEIGHT NUM_ERRORS NUM_ITERS DECODER [OPT]\n");
    fprintf(stderr, "e.g.:  ./mdpc-gf4 10 100 2293 37 88 200 0\n");
    fprintf(stderr, "\nParameter values:\n");
    fprintf(stderr, "NUM_KEYS:     positive integer, number of key pairs to generate\n");
    fprintf(stderr, "NUM_MSGS:     positive integer, number of messages to generate, encrypt and decrypt for every key pair\n");
    fprintf(stderr, "BLOCK_SIZE:   positive integer, size of cyclic matrix block, 2293 or 2339 is recommended\n");
    fprintf(stderr, "BLOCK_WEIGHT: positive integer, hamming weight of cyclic matrix block, 37 is recommended\n");
    fprintf(stderr, "NUM_ERRORS:   positive integer, number of errors in the error vector\n");
    fprintf(stderr, "NUM_ITERS:    positive integer, number of decoding iterations\n");
    fprintf(stderr, "DECODER:      positive integer, one of the options listed bellow\n");
    fprintf(stderr, "OPT:          optional positive integer, used as delta with decoders 2 and 3\n");

    fprintf(stderr, "\nPossible OPT settings:\n");
    fprintf(stderr, "With decoder 2, OPT as used as value of delta.\n");
    fprintf(stderr, "With decoder 3:\n");
    fprintf(stderr, "0 --> threshold function T0\n");
    fprintf(stderr, "1 --> threshold function T1\n");
    fprintf(stderr, "2 --> threshold function T2\n");
    fprintf(stderr, "3 --> threshold function T3\n");
    fprintf(stderr, "4 --> threshold function T4\n");
    fprintf(stderr, "5 --> threshold function T5\n");



    fprintf(stderr, "\nPossible decoders:\n");
    fprintf(stderr, "0 --> SF v1\n");
    fprintf(stderr, "1 --> SF v2  DEPRECATED! \n");
    fprintf(stderr, "2 --> SF with delta\n");
    fprintf(stderr, "3 --> SF with thr\n");
    fprintf(stderr, "4 --> SF BG\n");
}

int main(int argc, char ** argv) {
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
    size_t opt;
    num_keys = atol(argv[1]);
    num_messages = atol(argv[2]);
    block_size = atol(argv[3]);
    if (2293 != block_size && 2339 != block_size) {
        fprintf(stderr, "WARNING: recommended block size is 2293 or 2339! Provided value %zu is untested!\n", block_size);
    }
    block_weight = atol(argv[4]);
    if (37 != block_weight) {
        fprintf(stderr, "WARNING: recommended block weight is 37! Provided value %zu is untested!\n", block_weight);
    }
    num_errors = atol(argv[5]);
    if (num_errors < 84) {
        fprintf(stderr, "WARNING: recommended number of errors is 84! Provided value %zu is untested!\n", num_errors);
    }
    num_iterations = atol(argv[6]);
    decoder = atol(argv[7]);
    if (decoder > 4) {
        fprintf(stderr, "ERROR: possible decoders are 0-4! Provided value %zu is unsupported!\n", decoder);
        return -1;
    }
    if (decoder >= 2 && argc == 8) {
        fprintf(stderr, "ERROR: decoders 2, 3 and 4 require you to specify OPT! No OPT value was provided!\n");
        return -1;
    }
    if (1 == decoder) {
        fprintf(stderr, "WARNING: Decoder SF v2 is deprecated and should not be used!\n");
    }
    fprintf(stderr, "Params are: %zu %zu %zu %zu %zu %zu %zu", num_keys, num_messages, block_size, block_weight, num_errors, num_iterations, decoder);
    if (8 == argc) {
        fprintf(stderr, "\n");
        run_tests(stderr, num_keys, num_messages, block_size, block_weight, num_errors, num_iterations, decoder, 0);
    } else { // opt is provided
        opt = atol(argv[8]);
        if (3 == decoder && opt > 5) {
            fprintf(stderr, "ERROR: possible opt values for decoder 3 are 0-4! Provided value %zu is unsupported!\n", opt);
        }
        fprintf(stderr, " %zu\n", opt);
        run_tests(stderr, num_keys, num_messages, block_size, block_weight, num_errors, num_iterations, decoder, opt);
    }
    return 0;
}
#endif