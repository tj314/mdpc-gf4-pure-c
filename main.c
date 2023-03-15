#include <stdio.h>
#include "src/gf4.h"
#include "src/gf4_poly.h"
#include "src/contexts.h"
#include "src/enc.h"
#include "src/dec.h"


void test_load() {
    size_t block_size = 2339;
    size_t block_weight = 37;
    // size_t block_size = 7;
    // size_t block_weight = 3;
    encoding_context_t ec_gen, ec_load;
    decoding_context_t dc_gen, dc_load;
    fprintf(stderr, "Generating contexts!\n");
    contexts_init(&ec_gen, &dc_gen, block_size, block_weight);
    fprintf(stderr, "Contexts generated!\n");

    contexts_save("save.bin", &ec_gen, &dc_gen);
    contexts_load("save.bin", block_size, &ec_load, &dc_load);

    for (size_t i = 0; i < block_size; ++i) {
        if (ec_gen.second_block_G.coefficients[i] != ec_load.second_block_G.coefficients[i] ||
        dc_gen.h0.coefficients[i] != dc_load.h0.coefficients[i] ||
        dc_gen.h1.coefficients[i] != dc_load.h1.coefficients[i]) {
            fprintf(stderr, "coefficients do not match!\n");
            exit(-1);
        }
    }
    fprintf(stderr, "all ok!\n");
    contexts_deinit(&ec_gen, &dc_gen);
    contexts_deinit(&ec_load, &dc_load);
}

void test_inverse() {
    gf4_poly_t poly1 = gf4_poly_init_zero(10);
    gf4_poly_t poly2 = gf4_poly_init_zero(10);
    gf4_poly_t modulus = gf4_poly_init_zero(10);
    gf4_poly_t maybe_inverse = gf4_poly_init_zero(10);

    gf4_poly_set_coefficient(&poly1, 2, 1);
    gf4_poly_set_coefficient(&poly1, 1, 1);
    gf4_poly_set_coefficient(&poly1, 0, 1);

    gf4_poly_set_coefficient(&poly2, 4, 2);
    gf4_poly_set_coefficient(&poly2, 1, 2);

    gf4_poly_set_coefficient(&modulus, 0, 1);
    gf4_poly_set_coefficient(&modulus, 8, 1);

    gf4_poly_t * poly = &poly1;
    bool inverse_exists = gf4_poly_invert_slow(&maybe_inverse, poly, &modulus);


    printf("poly:    ");
    gf4_poly_pretty_print(poly, stderr, "\n");
    printf("modulus: ");
    gf4_poly_pretty_print(&modulus, stderr, "\n");

    if (inverse_exists) {
        printf("inverse: ");
        gf4_poly_pretty_print(&maybe_inverse, stderr, "\n");
        gf4_poly_t div = gf4_poly_init_zero(10);
        gf4_poly_t rem = gf4_poly_init_zero(10);
        gf4_poly_t m = gf4_poly_init_zero(20);
        gf4_poly_mul(&m, &maybe_inverse, poly);
        gf4_poly_div_rem(&div, &rem, &m, &modulus);
        printf("inv*poly %% modulus: ");
        gf4_poly_pretty_print(&rem, stderr, "\n");
        gf4_poly_deinit(&div);
        gf4_poly_deinit(&rem);
        gf4_poly_deinit(&m);
    } else {
        printf("inverse does not exist!\n");
    }

    gf4_poly_deinit(&poly1);
    gf4_poly_deinit(&poly2);
    gf4_poly_deinit(&modulus);
    gf4_poly_deinit(&maybe_inverse);
}

void test_encrypt_decrypt(size_t block_size, size_t block_weight, size_t num_errors, size_t num_iterations) {
    encoding_context_t ec;
    decoding_context_t dc;
    contexts_init(&ec, &dc, block_size, block_weight);

    gf4_poly_t msg = gf4_poly_init_zero(block_size);
    gf4_poly_t encrypted = gf4_poly_init_zero(2*block_size);
    gf4_poly_t decoded = gf4_poly_init_zero(2*block_size);

    random_weighted_gf4_poly(&msg, block_size, block_weight);
    enc_encrypt(&encrypted, &msg, num_errors, &ec);

    bool decoding_success = dec_decrypt(&decoded, &encrypted, dec_decode_symbol_flipping, num_iterations, &dc);
    if (decoding_success) {
        fprintf(stderr, "decryption success!\n");
        fprintf(stderr, "msg: ");
        gf4_poly_coeff_print(&msg, block_size, stderr, "\n");
        fprintf(stderr, "dec: ");
        gf4_poly_coeff_print(&decoded, block_size, stderr, "\n");
        bool matching = true;
        size_t num_misses = 0;
        for (size_t i = 0; i < block_size; ++i) {
            bool tmp = decoded.coefficients[i] == msg.coefficients[i];
            matching = matching && tmp;
            num_misses += (size_t)tmp;
        }
        if (matching) {
            fprintf(stderr, "decryption correct!\n");
        } else {
            fprintf(stderr, "decryption incorrect!!!\n");
            fprintf(stderr, "num misses: %zu\n", num_misses);
        }
    } else {
        fprintf(stderr, "decryption failure!\n");
    }

    gf4_poly_deinit(&msg);
    gf4_poly_deinit(&encrypted);
    gf4_poly_deinit(&decoded);
    contexts_deinit(&ec, &dc);
}

void test_encrypt_decrypt_from_file(const char * filepath,
                                    size_t block_size, size_t block_weight, size_t num_errors, size_t num_iterations) {
    encoding_context_t ec;
    decoding_context_t dc;
    contexts_load(filepath, block_size, &ec, &dc);

    gf4_poly_t msg = gf4_poly_init_zero(block_size);
    gf4_poly_t encrypted = gf4_poly_init_zero(2*block_size);
    gf4_poly_t decoded = gf4_poly_init_zero(2*block_size);

    random_weighted_gf4_poly(&msg, block_size, block_weight);
    enc_encrypt(&encrypted, &msg, num_errors, &ec);

    bool decoding_success = dec_decrypt(&decoded, &encrypted, dec_decode_symbol_flipping, num_iterations, &dc);
    if (decoding_success) {
        fprintf(stderr, "decryption success!\n");
        fprintf(stderr, "msg: ");
        gf4_poly_coeff_print(&msg, block_size, stderr, "\n");
        fprintf(stderr, "dec: ");
        gf4_poly_coeff_print(&decoded, block_size, stderr, "\n");
        bool matching = true;
        size_t num_misses = 0;
        for (size_t i = 0; i < block_size; ++i) {
            bool tmp = decoded.coefficients[i] == msg.coefficients[i];
            matching = matching && tmp;
            num_misses += (size_t)tmp;
        }
        if (matching) {
            fprintf(stderr, "decryption correct!\n");
        } else {
            fprintf(stderr, "decryption incorrect!!!\n");
            fprintf(stderr, "num misses: %zu\n", num_misses);
        }
    } else {
        fprintf(stderr, "decryption failure!\n");
    }

    gf4_poly_deinit(&msg);
    gf4_poly_deinit(&encrypted);
    gf4_poly_deinit(&decoded);
    contexts_deinit(&ec, &dc);
}

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



int main(void) {
    // test_inverse();
    // test_load();

    size_t block_size = 2339;
    size_t block_weight = 37;
    size_t num_errors = 84;
    size_t num_iterations = 150;
    // size_t block_size = 7;
    // size_t block_weight = 3;
    // size_t num_errors = 1;
    // size_t num_iterations = 10;
    // test_encrypt_decrypt(block_size, block_weight, num_errors, num_iterations);

    gen_keys(block_size, block_weight);

    return 0;
}
