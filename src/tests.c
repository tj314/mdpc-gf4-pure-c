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

void assert_compare_coeffs(gf4_t * coeffs1, gf4_t * coeffs2, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        assert(coeffs1[i] == coeffs2[i]);
    }
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

    // test 1
    gf4_vector_t vector = gf4_vector_init_with_length(10, 10);
    assert(0 == utils_hamming_weight(&vector));
    assert(10 == vector.capacity);

    // test 2
    for (size_t i = 0; i < 4; ++i) {
        size_t deg = random_from_range(0, 9);
        while (0 != vector.array[deg]) { // this index is already set, choose another one
            deg = random_from_range(0, 9);
        }
        gf4_t val = random_from_range(1, GF4_MAX_VALUE);
        vector.array[deg] = val;
    }
    assert(4 == utils_hamming_weight(&vector));
    assert(10 == vector.capacity);

    // cleanup
    gf4_vector_deinit(&vector);
    print_OK();
}

// gf4_poly
void test_gf4_poly_init_zero(){
    fprintf(stderr, "test_gf4_poly_init_zero: FAIL\n");
    // print_OK();
}

void test_gf4_poly_zero_out(){
    fprintf(stderr, "test_gf4_poly_zero_out: FAIL\n");
    // print_OK();
}

void test_gf4_poly_deinit(){
    fprintf(stderr, "test_gf4_poly_deinit: FAIL\n");
    // print_OK();
}

void test_gf4_poly_get_coefficient(){
    fprintf(stderr, "test_gf4_poly_get_coefficient: FAIL\n");
    // print_OK();
}

void test_gf4_poly_set_coefficient(){
    fprintf(stderr, "test_gf4_poly_set_coefficient: FAIL\n");
    // print_OK();
}

void test_gf4_poly_add(){
    fprintf(stderr, "test_gf4_poly_add: FAIL\n");
    // print_OK();
}

void test_gf4_poly_add_inplace(){
    fprintf(stderr, "test_gf4_poly_add_inplace: FAIL\n");
    // print_OK();
}

void test_gf4_poly_add_ax_to_deg_inplace(){
    fprintf(stderr, "test_gf4_poly_add_ax_to_deg_inplace: FAIL\n");
    // print_OK();
}

void test_gf4_poly_mul(){
    fprintf(stderr, "test_gf4_poly_mul: FAIL\n");
    // print_OK();
}

void test_gf4_poly_div_x_to_deg(){
    fprintf(stderr, "test_gf4_poly_div_x_to_deg: FAIL\n");
    // print_OK();
}

void test_gf4_poly_div_x_to_deg_inplace(){
    fprintf(stderr, "test_gf4_poly_div_x_to_deg_inplace: FAIL\n");
    // print_OK();
}

void test_gf4_poly_div_rem(){
    fprintf(stderr, "test_gf4_poly_div_rem: FAIL\n");
    // print_OK();
}

void test_gf4_poly_invert_slow(){
    fprintf(stderr, "test_gf4_poly_invert_slow: FAIL\n");
    // print_OK();
}

void test_gf4_poly_is_zero(){
    fprintf(stderr, "test_gf4_poly_is_zero: FAIL\n");
    // print_OK();
}

void test_gf4_poly_equal(){
    fprintf(stderr, "test_gf4_poly_equal: FAIL\n");
    // print_OK();
}

void test_gf4_poly_cyclic_shift_right_inplace() {
    fprintf(stderr, "test_gf4_poly_cyclic_shift_right_inplace: ");
    // setup
    gf4_poly_t poly = gf4_poly_init_zero(5);
    gf4_poly_set_coefficient(&poly, 0, 1);
    gf4_poly_set_coefficient(&poly, 1, 1);
    gf4_poly_set_coefficient(&poly, 3, 1);
    gf4_t expected_poly[5] = {1,1,0,1,0};
    assert(5 == poly.capacity);
    assert(3 == gf4_poly_get_degree(&poly));
    assert_compare_coeffs(poly.array, expected_poly, 5);

    // test 1
    // cyclic shift of 1 + x + x^3 with size=4 and capacity=5: [1, 1, 0, 1, 0] -> [1, 1, 1, 0, 0] -> 1 + x + x^2
    gf4_poly_t poly1 = gf4_poly_clone(&poly);
    gf4_t expected_poly1_after_shift[5] = {1,1,1,0,0};
    assert(5 == poly1.capacity);
    assert(3 == gf4_poly_get_degree(&poly1));
    assert_compare_coeffs(poly1.array, expected_poly, 5);

    gf4_poly_cyclic_shift_right_inplace(&poly1, 4);
    assert(5 == poly1.capacity);
    assert(3 == gf4_poly_get_degree(&poly1));
    assert(5 == poly.capacity);
    assert(3 == gf4_poly_get_degree(&poly));
    assert_compare_coeffs(poly.array, expected_poly, 5);
    assert_compare_coeffs(poly1.array, expected_poly1_after_shift, 5);
    gf4_poly_deinit(&poly1);

    // test 2
    // cyclic shift of 1 + x + x^3 with size=5 and capacity=5: [1, 1, 0, 1, 0] -> [0, 1, 1, 0, 1] -> x + x^2 + x^4
    gf4_poly_t poly2 = gf4_poly_clone(&poly);
    assert(5 == poly2.capacity);
    assert(3 == gf4_poly_get_degree(&poly2));
    assert_compare_coeffs(poly2.array, expected_poly, 5);

    gf4_t expected_poly2_after_shift[5] = {0, 1, 1, 0, 1};
    gf4_poly_cyclic_shift_right_inplace(&poly2, 5);
    assert(5 == poly2.capacity);
    assert(3 == gf4_poly_get_degree(&poly2));
    assert_compare_coeffs(poly2.array, expected_poly2_after_shift, 5);
    assert(5 == poly.capacity);
    assert(3 == gf4_poly_get_degree(&poly));
    assert_compare_coeffs(poly.array, expected_poly, 5);
    gf4_poly_deinit(&poly2);

    // cleanup
    gf4_poly_deinit(&poly);
    print_OK();
}

void test_gf4_poly_get_degree() {
    fprintf(stderr, "test_gf4_poly_get_degree: FAIL\n");
    // print_OK();
}
void test_gf4_poly_adjust_degree() {
    fprintf(stderr, "test_gf4_poly_adjust_degree: FAIL\n");
    // print_OK();
}
void test_gf4_poly_pretty_print() {
    fprintf(stderr, "test_gf4_poly_pretty_print: FAIL\n");
    // print_OK();
}
void test_gf4_poly_clone() {
    fprintf(stderr, "test_gf4_poly_clone: FAIL\n");
    // print_OK();
}
void test_gf4_poly_copy() {
    fprintf(stderr, "test_gf4_poly_copy: FAIL\n");
    // print_OK();
}


// gf4_matrix
void test_gf4_square_matrix_init_cyclic_matrix() {
    // setup
    fprintf(stderr, "test_gf4_square_matrix_init_cyclic_matrix: ");
    gf4_poly_t first_row = gf4_poly_init_zero(5);
    gf4_poly_set_coefficient(&first_row, 0, 1);
    gf4_poly_set_coefficient(&first_row, 1, 2);
    gf4_poly_set_coefficient(&first_row, 4, 3);

    gf4_t expected_matrix[5][5] = {
            {1, 2, 0, 0, 3},
            {3, 1, 2, 0, 0},
            {0, 3, 1, 2, 0},
            {0, 0, 3, 1, 2},
            {2, 0, 0, 3, 1}
    };

    gf4_matrix_t matrix = gf4_matrix_init_cyclic_matrix(&first_row, 5);
    assert(5 == matrix.num_rows);
    for (size_t i = 0; i < 5; ++i) {
        assert_compare_coeffs(matrix.rows[i], expected_matrix[i], 5);
    }

    // cleanup
    gf4_poly_deinit(&first_row);
    gf4_matrix_deinit(&matrix);
    print_OK();
}

void test_gf4_matrix_gaussian_elimination_inplace() {
    fprintf(stderr, "test_gf4_matrix_gaussian_elimination_inplace: ");
    {
        gf4_poly_t first_row = gf4_poly_init_zero(4);
        gf4_poly_set_coefficient(&first_row, 0, 1);
        gf4_poly_set_coefficient(&first_row, 1, 2);
        gf4_poly_set_coefficient(&first_row, 2, 1);
        gf4_poly_set_coefficient(&first_row, 3, 3);
        gf4_matrix_t matrix = gf4_matrix_init_cyclic_matrix(&first_row, 4);
        gf4_matrix_gaussian_elimination_inplace(&matrix);
        gf4_t expected_matrix[4][4] = {
                {2, 0, 0, 0},
                {0, 3, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 3}
        };
        // gf4_matrix_pretty_print(&matrix);
        assert(4 == matrix.num_rows);
        assert(4 == matrix.num_cols);
        for (size_t i = 0; i < 4; ++i) {
            assert_compare_coeffs(matrix.rows[i], expected_matrix[i], 4);
        }
        gf4_poly_deinit(&first_row);
    }
    {
        gf4_poly_t first_row = gf4_poly_init_zero(4);
        gf4_poly_set_coefficient(&first_row, 0, 1);
        gf4_poly_set_coefficient(&first_row, 2, 1);
        gf4_matrix_t matrix = gf4_matrix_init_cyclic_matrix(&first_row, 4);
        gf4_matrix_gaussian_elimination_inplace(&matrix);
        gf4_t expected_matrix[4][4] = {
                {1, 0, 1, 0},
                {0, 1, 0, 1},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        };
        // gf4_matrix_pretty_print(&matrix);
        assert(4 == matrix.num_rows);
        assert(4 == matrix.num_cols);
        for (size_t i = 0; i < 4; ++i) {
            assert_compare_coeffs(matrix.rows[i], expected_matrix[i], 4);
        }
        gf4_poly_deinit(&first_row);
    }
    print_OK();
}

void test_gf4_matrix_solve_homogenous_linear_system() {
    fprintf(stderr, "test_gf4_matrix_solve_homogenous_linear_system: ");

    // test 1
    {
        gf4_matrix_t matrix;
        matrix.num_rows = 1;
        matrix.num_cols = 3;
        matrix.rows = malloc(sizeof(gf4_t *));
        assert(NULL != matrix.rows);
        matrix.rows[0] = malloc(3 * sizeof(gf4_t));
        assert(NULL != matrix.rows[0]);
        matrix.rows[0][0] = 2;
        matrix.rows[0][1] = 1;
        matrix.rows[0][2] = 0;

        gf4_matrix_gaussian_elimination_inplace(&matrix);
        assert(1 == matrix.num_rows);
        assert(3 == matrix.num_cols);
        assert(2 == matrix.rows[0][0]);
        assert(1 == matrix.rows[0][1]);
        assert(0 == matrix.rows[0][2]);

        gf4_matrix_t solutions = gf4_matrix_solve_homogenous_linear_system(&matrix);
        gf4_t expected_solutions[16][3] = {
                {0, 0, 0},
                {0, 0, 1},
                {0, 0, 2},
                {0, 0, 3},
                {3, 1, 0},
                {3, 1, 1},
                {3, 1, 2},
                {3, 1, 3},
                {1, 2, 0},
                {1, 2, 1},
                {1, 2, 2},
                {1, 2, 3},
                {2, 3, 0},
                {2, 3, 1},
                {2, 3, 2},
                {2, 3, 3}
        };
        assert(16 == solutions.num_rows);
        assert(3 == solutions.num_cols);
        for (size_t i = 0; i < 16; ++i) {
            assert_compare_coeffs(solutions.rows[i], expected_solutions[i], 3);
        }

        gf4_matrix_deinit(&matrix);
        gf4_matrix_deinit(&solutions);
    }

    // test 2
    {
        gf4_poly_t first_row = gf4_poly_init_zero(3);
        gf4_poly_set_coefficient(&first_row, 0, 2);
        gf4_matrix_t equations = gf4_matrix_init_cyclic_matrix(&first_row, 3);

        gf4_t expected_equations[3][3] = {
                {2, 0, 0},
                {0, 2, 0},
                {0, 0, 2}
        };

        assert(3 == equations.num_cols);
        assert(3 == equations.num_rows);
        for (size_t i = 0; i < 3; ++i) {
            assert_compare_coeffs(equations.rows[i], expected_equations[i], 3);
        }

        gf4_matrix_gaussian_elimination_inplace(&equations);
        assert(3 == equations.num_cols);
        assert(3 == equations.num_rows);
        for (size_t i = 0; i < 3; ++i) {
            assert_compare_coeffs(equations.rows[i], expected_equations[i], 3);
        }

        gf4_matrix_t solutions = gf4_matrix_solve_homogenous_linear_system(&equations);
        gf4_t expected_solution[3] = {0, 0, 0};
        assert(1 == solutions.num_rows);
        assert(3 == solutions.num_cols);
        assert_compare_coeffs(solutions.rows[0], expected_solution, 3);

        gf4_matrix_deinit(&solutions);
        gf4_matrix_deinit(&equations);
        gf4_poly_deinit(&first_row);
    }

    // test 3
    {
        gf4_matrix_t matrix;
        matrix.num_rows = 1;
        matrix.num_cols = 3;
        matrix.rows = malloc(sizeof(gf4_t *));
        assert(NULL != matrix.rows);
        matrix.rows[0] = malloc(3 * sizeof(gf4_t));
        assert(NULL != matrix.rows[0]);
        matrix.rows[0][0] = 0;
        matrix.rows[0][1] = 1;
        matrix.rows[0][2] = 0;

        gf4_matrix_gaussian_elimination_inplace(&matrix);
        assert(1 == matrix.num_rows);
        assert(3 == matrix.num_cols);
        assert(0 == matrix.rows[0][0]);
        assert(1 == matrix.rows[0][1]);
        assert(0 == matrix.rows[0][2]);

        gf4_matrix_t solutions = gf4_matrix_solve_homogenous_linear_system(&matrix);
        gf4_t expected_solutions[16][3] = {
                {0, 0, 0},
                {0, 0, 1},
                {0, 0, 2},
                {0, 0, 3},
                {1, 0, 0},
                {1, 0, 1},
                {1, 0, 2},
                {1, 0, 3},
                {2, 0, 0},
                {2, 0, 1},
                {2, 0, 2},
                {2, 0, 3},
                {3, 0, 0},
                {3, 0, 1},
                {3, 0, 2},
                {3, 0, 3}
        };
        assert(16 == solutions.num_rows);
        assert(3 == solutions.num_cols);
        for (size_t i = 0; i < 16; ++i) {
            assert_compare_coeffs(solutions.rows[i], expected_solutions[i], 3);
        }

        gf4_matrix_deinit(&matrix);
        gf4_matrix_deinit(&solutions);
    }

    // test 4
    {
        gf4_matrix_t matrix;
        matrix.num_rows = 4;
        matrix.num_cols = 4;
        matrix.rows = malloc(4* sizeof(gf4_t *));
        assert(NULL != matrix.rows);
        for (size_t i = 0; i < 4; ++i) {
            matrix.rows[i] = calloc(4, sizeof(gf4_t));
            assert(NULL != matrix.rows[0]);
        }
        matrix.rows[0][0] = 1;
        matrix.rows[0][1] = 2;
        matrix.rows[1][2] = 1;
        matrix.rows[2][3] = 3;

        gf4_matrix_gaussian_elimination_inplace(&matrix);
        assert(4 == matrix.num_rows);
        assert(4 == matrix.num_cols);
        assert(1 == matrix.rows[0][0]);
        assert(2 == matrix.rows[0][1]);
        assert(1 == matrix.rows[1][2]);
        assert(3 == matrix.rows[2][3]);

        gf4_matrix_t solutions = gf4_matrix_solve_homogenous_linear_system(&matrix);
        gf4_t expected_solutions[4][4] = {
                {0, 0, 0, 0},
                {2, 1, 0, 0},
                {3, 2, 0, 0},
                {1, 3, 0, 0}
        };
        assert(4 == solutions.num_rows);
        assert(4 == solutions.num_cols);
        for (size_t i = 0; i < 4; ++i) {
            assert_compare_coeffs(solutions.rows[i], expected_solutions[i], 4);
        }

        gf4_matrix_deinit(&matrix);
        gf4_matrix_deinit(&solutions);
    }

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
    sprintf(filename, "test-file-%lu.txt", (unsigned long)time(NULL));

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
    remove(filename);
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
    assert(2 == ec.second_block_G.array[0]);
    assert(1 == ec.second_block_G.array[1]);
    assert(3 == ec.second_block_G.array[2]);
    assert(3 == ec.block_size);
    assert(3 == ec.second_block_G.capacity);
    assert(6 == encoded.capacity);
    for (size_t i = 0; i < 6; ++i) assert(0 == encoded.array[i]);

    // test 2
    gf4_poly_set_coefficient(&msg, 0, 1);
    gf4_poly_set_coefficient(&msg, 2, 2);
    enc_encode(&encoded, &msg, &ec);
    assert(2 == ec.second_block_G.array[0]);
    assert(1 == ec.second_block_G.array[1]);
    assert(3 == ec.second_block_G.array[2]);
    assert(3 == ec.block_size);
    assert(3 == ec.second_block_G.capacity);
    assert(3 == msg.capacity);
    assert(1 == msg.array[0]);
    assert(0 == msg.array[1]);
    assert(2 == msg.array[2]);
    assert(0 == memcmp(encoded.array, msg.array, 3));
    assert(3 == encoded.array[3]);
    assert(1 == encoded.array[4]);
    assert(2 == encoded.array[5]);

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
        random_gf4_vector(&msg, 3);
        enc_encode(&encoded, &msg, &ec);
        enc_encrypt(&encrypted, &msg, 1, &ec);

        size_t differences = 0;
        for (size_t j = 0; j < 6; ++j) {
            differences += (encoded.array[j] != encrypted.array[j]);
        }

        assert(1 == differences);

        assert(0 == ec.second_block_G.array[0]);
        assert(2 == ec.second_block_G.array[1]);
        assert(3 == ec.second_block_G.array[2]);
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
    assert(0 == syndrome.array[0]);
    assert(0 == syndrome.array[1]);
    assert(0 == syndrome.array[2]);
    assert(0 == vec.array[0]);
    assert(0 == vec.array[1]);
    assert(0 == vec.array[2]);
    assert(0 == vec.array[3]);
    assert(0 == vec.array[4]);
    assert(0 == vec.array[5]);
    assert(1 == dc.h0.array[0]);
    assert(0 == dc.h0.array[1]);
    assert(0 == dc.h0.array[2]);
    assert(2 == dc.h1.array[0]);
    assert(0 == dc.h1.array[1]);
    assert(3 == dc.h1.array[2]);
    assert(3 == dc.block_size);
    assert(0 == gf4_poly_get_degree(&dc.h0));
    assert(2 == gf4_poly_get_degree(&dc.h1));

    // test 2
    gf4_poly_set_coefficient(&vec, 0, 1);
    gf4_poly_set_coefficient(&vec, 1, 2);
    gf4_poly_set_coefficient(&vec, 4, 1);
    gf4_poly_set_coefficient(&vec, 5, 3);

    dec_calculate_syndrome(&syndrome, &vec, &dc);

    assert(3 == syndrome.array[0]);
    assert(0 == syndrome.array[1]);
    assert(2 == syndrome.array[2]);
    assert(1 == vec.array[0]);
    assert(2 == vec.array[1]);
    assert(0 == vec.array[2]);
    assert(0 == vec.array[3]);
    assert(1 == vec.array[4]);
    assert(3 == vec.array[5]);
    assert(1 == dc.h0.array[0]);
    assert(0 == dc.h0.array[1]);
    assert(0 == dc.h0.array[2]);
    assert(2 == dc.h1.array[0]);
    assert(0 == dc.h1.array[1]);
    assert(3 == dc.h1.array[2]);
    assert(3 == dc.block_size);
    assert(0 == gf4_poly_get_degree(&dc.h0));
    assert(2 == gf4_poly_get_degree(&dc.h1));

    // cleanup
    gf4_poly_deinit(&vec);
    gf4_poly_deinit(&dc.h0);
    gf4_poly_deinit(&dc.h1);
    gf4_poly_deinit(&syndrome);
    print_OK();
}