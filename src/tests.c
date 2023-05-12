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

#include "tests.h"

void print_OK() {
    fprintf(stderr, "\x1B[32m");
    fprintf(stderr, "OK\n");
    fprintf(stderr, "\x1B[0m");
}


// gf4
void test_gf4_is_in_range() {
    fprintf(stderr, "test_gf4_is_in_range: ");
    assert(true == gf4_is_in_range(0));
    assert(true == gf4_is_in_range(1));
    assert(true == gf4_is_in_range(2));
    assert(true == gf4_is_in_range(3));
    assert(false == gf4_is_in_range(-1));
    assert(false == gf4_is_in_range(4));
    print_OK();
}

void test_gf4_to_str() {
    fprintf(stderr, "test_gf4_to_str: ");
    assert(0 == strcmp("0", gf4_to_str(0)));
    assert(0 == strcmp("1", gf4_to_str(1)));
    assert(0 == strcmp("a", gf4_to_str(2)));
    assert(0 == strcmp("(a+1)", gf4_to_str(3)));
    print_OK();
}

void test_gf4_add() {
    fprintf(stderr, "test_gf4_add: ");
    assert(0 == gf4_add(0, 0));
    assert(1 == gf4_add(1, 0));
    assert(2 == gf4_add(2, 0));
    assert(3 == gf4_add(3, 0));
    assert(1 == gf4_add(0, 1));
    assert(0 == gf4_add(1, 1));
    assert(3 == gf4_add(2, 1));
    assert(2 == gf4_add(3, 1));
    assert(2 == gf4_add(0, 2));
    assert(3 == gf4_add(1, 2));
    assert(0 == gf4_add(2, 2));
    assert(1 == gf4_add(3, 2));
    assert(3 == gf4_add(0, 3));
    assert(2 == gf4_add(1, 3));
    assert(1 == gf4_add(2, 3));
    assert(0 == gf4_add(3, 3));
    print_OK();
}

void test_gf4_mul() {
    fprintf(stderr, "test_gf4_mul: ");
    assert(0 == gf4_mul(0, 0));
    assert(0 == gf4_mul(1, 0));
    assert(0 == gf4_mul(2, 0));
    assert(0 == gf4_mul(3, 0));
    assert(0 == gf4_mul(0, 1));
    assert(1 == gf4_mul(1, 1));
    assert(2 == gf4_mul(2, 1));
    assert(3 == gf4_mul(3, 1));
    assert(0 == gf4_mul(0, 2));
    assert(2 == gf4_mul(1, 2));
    assert(3 == gf4_mul(2, 2));
    assert(1 == gf4_mul(3, 2));
    assert(0 == gf4_mul(0, 3));
    assert(3 == gf4_mul(1, 3));
    assert(1 == gf4_mul(2, 3));
    assert(2 == gf4_mul(3, 3));
    print_OK();
}

void test_gf4_div() {
    fprintf(stderr, "test_gf4_div: ");
    assert(0 == gf4_div(0, 1));
    assert(1 == gf4_div(1, 1));
    assert(2 == gf4_div(2, 1));
    assert(3 == gf4_div(3, 1));
    assert(0 == gf4_div(0, 2));
    assert(3 == gf4_div(1, 2));
    assert(1 == gf4_div(2, 2));
    assert(2 == gf4_div(3, 2));
    assert(0 == gf4_div(0, 3));
    assert(2 == gf4_div(1, 3));
    assert(3 == gf4_div(2, 3));
    assert(1 == gf4_div(3, 3));
    print_OK();
}

// utils
void test_utils_hamming_weight() {
    fprintf(stderr, "test_utils_hamming_weight: ");
    gf4_poly_t p = gf4_poly_init_zero(10);
    assert(0 == utils_hamming_weight(&p));
    assert(10 == p.capacity);
    for (size_t i = 0; i < 4; ++i) {
        size_t deg = random_from_range(0, 9);
        while (0 != p.coefficients[deg]) {
            deg = random_from_range(0, 9);
        }
        gf4_t val = random_from_range(1, GF4_MAX_VALUE);
        gf4_poly_set_coefficient(&p, deg, val);
    }
    assert(4 == utils_hamming_weight(&p));
    assert(10 == p.capacity);
    gf4_poly_deinit(&p);
    print_OK();
}

// contexts
void test_contexts_init() {
    fprintf(stderr, "test_contexts_init: ");
    // setup
    const size_t block_size = 2339;
    const size_t capacity = block_size + 1;
    gf4_poly_t inv = gf4_poly_init_zero(2*capacity);
    gf4_poly_t tmp = gf4_poly_init_zero(2*capacity);
    gf4_poly_t div = gf4_poly_init_zero(2*capacity);
    gf4_poly_t rem = gf4_poly_init_zero(2*capacity);
    gf4_poly_t modulus = gf4_poly_init_zero(capacity);
    gf4_poly_set_coefficient(&modulus, 0, 1);
    gf4_poly_set_coefficient(&modulus, block_size, 1);

    // tests
    for (size_t i = 0; i < 10; ++i) {
        encoding_context_t ec;
        decoding_context_t dc;
        contexts_init(&ec, &dc, block_size, 37);

        // Test1: h1 must be invertible
        assert(gf4_poly_invert_slow(&inv, &dc.h1, &modulus));

        // Test2:
        // H*G^T must be 0
        // H = (H0|H1), G=(I|-(H0*H1^-1)^T)
        // H*G^T = H0 - H1*H0*H1^-1 ---> dc.h0 - dc.h1*ec.second_block_G
        gf4_poly_mul(&tmp, &dc.h1, &ec.second_block_G);
        gf4_poly_div_rem(&div, &rem, &tmp, &modulus);
        gf4_poly_add_inplace(&rem, &dc.h0);
        assert(gf4_poly_is_zero(&rem));

        gf4_poly_zero_out(&div);
        gf4_poly_zero_out(&rem);
        gf4_poly_zero_out(&tmp);
        gf4_poly_zero_out(&inv);
        contexts_deinit(&ec, &dc);
    }

    // cleanup
    gf4_poly_deinit(&rem);
    gf4_poly_deinit(&div);
    gf4_poly_deinit(&inv);
    gf4_poly_deinit(&modulus);
    gf4_poly_deinit(&tmp);
    print_OK();
}

void test_contexts_save_load() {
    fprintf(stderr, "test_contexts_save_load: ");
    // setup
    const size_t block_size = 2339;
    const size_t block_weight = 37;
    encoding_context_t ec_gen, ec_load;
    decoding_context_t dc_gen, dc_load;
    contexts_init(&ec_gen, &dc_gen, block_size, block_weight);
    char filename[100] = {0};
    sprintf(filename, "test-file-%lu.bin", (unsigned long)time(NULL));

    // test: save, load, compare to the original
    contexts_save(filename, &ec_gen, &dc_gen);
    contexts_load(filename, block_size, &ec_load, &dc_load);
    assert(ec_gen.block_size == ec_load.block_size);
    assert(dc_gen.block_size == dc_load.block_size);
    assert(gf4_poly_equal(&ec_gen.second_block_G, &ec_load.second_block_G));
    assert(gf4_poly_equal(&dc_gen.h0, &dc_load.h0));
    assert(gf4_poly_equal(&dc_gen.h1, &dc_load.h1));

    // cleanup
    contexts_deinit(&ec_gen, &dc_gen);
    contexts_deinit(&ec_load, &dc_load);
    print_OK();
}

// enc
void test_enc_encode() {
    fprintf(stderr, "test_enc_encode: ");
    // setup
    encoding_context_t ec;
    ec.block_size = 3;
    ec.second_block_G = gf4_poly_init_zero(ec.block_size);
    gf4_poly_t msg = gf4_poly_init_zero(ec.block_size);
    gf4_poly_t encoded = gf4_poly_init_zero(2*ec.block_size);
    gf4_poly_set_coefficient(&ec.second_block_G, 0, 2);
    gf4_poly_set_coefficient(&ec.second_block_G, 1, 1);
    gf4_poly_set_coefficient(&ec.second_block_G, 2, 3);

    // test 1
    enc_encode(&encoded, &msg, &ec);
    assert(2 == ec.second_block_G.coefficients[0]);
    assert(1 == ec.second_block_G.coefficients[1]);
    assert(3 == ec.second_block_G.coefficients[2]);
    assert(3 == ec.block_size);
    assert(3 == ec.second_block_G.capacity);
    assert(6 == encoded.capacity);
    for (size_t i = 0; i < 6; ++i) assert(0 == encoded.coefficients[i]);

    // test 2
    gf4_poly_set_coefficient(&msg, 0, 1);
    gf4_poly_set_coefficient(&msg, 2, 2);
    enc_encode(&encoded, &msg, &ec);
    assert(2 == ec.second_block_G.coefficients[0]);
    assert(1 == ec.second_block_G.coefficients[1]);
    assert(3 == ec.second_block_G.coefficients[2]);
    assert(3 == ec.block_size);
    assert(3 == ec.second_block_G.capacity);
    assert(3 == msg.capacity);
    assert(1 == msg.coefficients[0]);
    assert(0 == msg.coefficients[1]);
    assert(2 == msg.coefficients[2]);
    assert(0 == memcmp(encoded.coefficients, msg.coefficients, 3));
    assert(3 == encoded.coefficients[3]);
    assert(1 == encoded.coefficients[4]);
    assert(2 == encoded.coefficients[5]);

    // cleanup
    gf4_poly_deinit(&ec.second_block_G);
    gf4_poly_deinit(&msg);
    gf4_poly_deinit(&encoded);
    print_OK();
}

void test_enc_encrypt() {
    fprintf(stderr, "test_enc_encrypt: ");
    // setup
    encoding_context_t ec;
    ec.block_size = 3;
    ec.second_block_G = gf4_poly_init_zero(ec.block_size);
    gf4_poly_t msg = gf4_poly_init_zero(ec.block_size);
    gf4_poly_t encoded = gf4_poly_init_zero(2*ec.block_size);
    gf4_poly_t encrypted = gf4_poly_init_zero(2*ec.block_size);
    gf4_poly_set_coefficient(&ec.second_block_G, 1, 2);
    gf4_poly_set_coefficient(&ec.second_block_G, 2, 3);

    // 5 runs of test
    for (size_t i = 0; i < 5; ++i) {
        random_gf4_poly(&msg, 3);
        enc_encode(&encoded, &msg, &ec);
        enc_encrypt(&encrypted, &msg, 1, &ec);

        size_t differences = 0;
        for (size_t j = 0; j < 6; ++j) {
            differences += (encoded.coefficients[j] != encrypted.coefficients[j]);
        }

        assert(1 == differences);

        assert(0 == ec.second_block_G.coefficients[0]);
        assert(2 == ec.second_block_G.coefficients[1]);
        assert(3 == ec.second_block_G.coefficients[2]);
        assert(3 == ec.block_size);
        assert(3 == ec.second_block_G.capacity);
        assert(6 == encoded.capacity);
        assert(6 == encrypted.capacity);
    }
    gf4_poly_deinit(&ec.second_block_G);
    gf4_poly_deinit(&msg);
    gf4_poly_deinit(&encoded);
    gf4_poly_deinit(&encrypted);
    print_OK();
}


// dec
void test_dec_calculate_syndrome() {
    fprintf(stderr, "test_dec_calculate_syndrome: ");
    // setup
    decoding_context_t dc;
    dc.block_size = 3;
    dc.h0 = gf4_poly_init_zero(dc.block_size);
    dc.h1 = gf4_poly_init_zero(dc.block_size);
    gf4_poly_set_coefficient(&dc.h0, 0, 1);
    gf4_poly_set_coefficient(&dc.h1, 0, 2);
    gf4_poly_set_coefficient(&dc.h1, 2, 3);
    gf4_poly_t syndrome = gf4_poly_init_zero(dc.block_size);

    // test 1
    gf4_poly_t vec = gf4_poly_init_zero(2*dc.block_size);
    dec_calculate_syndrome(&syndrome, &vec, &dc);
    assert(0 == syndrome.coefficients[0]);
    assert(0 == syndrome.coefficients[1]);
    assert(0 == syndrome.coefficients[2]);
    assert(0 == vec.coefficients[0]);
    assert(0 == vec.coefficients[1]);
    assert(0 == vec.coefficients[2]);
    assert(0 == vec.coefficients[3]);
    assert(0 == vec.coefficients[4]);
    assert(0 == vec.coefficients[5]);
    assert(1 == dc.h0.coefficients[0]);
    assert(0 == dc.h0.coefficients[1]);
    assert(0 == dc.h0.coefficients[2]);
    assert(2 == dc.h1.coefficients[0]);
    assert(0 == dc.h1.coefficients[1]);
    assert(3 == dc.h1.coefficients[2]);
    assert(3 == dc.block_size);
    assert(0 == dc.h0.degree);
    assert(2 == dc.h1.degree);

    // test 2
    gf4_poly_set_coefficient(&vec, 0, 1);
    gf4_poly_set_coefficient(&vec, 1, 2);
    gf4_poly_set_coefficient(&vec, 4, 1);
    gf4_poly_set_coefficient(&vec, 5, 3);

    dec_calculate_syndrome(&syndrome, &vec, &dc);

    assert(3 == syndrome.coefficients[0]);
    assert(0 == syndrome.coefficients[1]);
    assert(2 == syndrome.coefficients[2]);
    assert(1 == vec.coefficients[0]);
    assert(2 == vec.coefficients[1]);
    assert(0 == vec.coefficients[2]);
    assert(0 == vec.coefficients[3]);
    assert(1 == vec.coefficients[4]);
    assert(3 == vec.coefficients[5]);
    assert(1 == dc.h0.coefficients[0]);
    assert(0 == dc.h0.coefficients[1]);
    assert(0 == dc.h0.coefficients[2]);
    assert(2 == dc.h1.coefficients[0]);
    assert(0 == dc.h1.coefficients[1]);
    assert(3 == dc.h1.coefficients[2]);
    assert(3 == dc.block_size);
    assert(0 == dc.h0.degree);
    assert(2 == dc.h1.degree);

    // cleanup
    gf4_poly_deinit(&vec);
    gf4_poly_deinit(&dc.h0);
    gf4_poly_deinit(&dc.h1);
    gf4_poly_deinit(&syndrome);
    print_OK();
}