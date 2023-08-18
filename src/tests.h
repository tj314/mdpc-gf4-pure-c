/**
 *  @file   tests.h
 *  @brief  Unit tests.
 *  @author Tomáš Vavro
 *  @date   2023-05-12
 ***********************************************/

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

#ifndef MDPC_GF4_TESTS_H
#define MDPC_GF4_TESTS_H

#include <stdio.h>
#include <stdlib.h>
#include "dec.h"
#include "contexts.h"
#include "enc.h"
#include "gf4.h"
#include "gf4_array.h"
#include "gf4_poly.h"
#include "gf4_matrix.h"
#include "random.h"
#include "utils.h"

// TESTS
// gf4
void test_gf4_is_in_range();
void test_gf4_to_str();
void test_gf4_add();
void test_gf4_mul();
void test_gf4_div();

// array
void test_gf4_array_hamming_weight();

// poly
void test_gf4_poly_init_zero();
void test_gf4_poly_zero_out();
void test_gf4_poly_deinit();
void test_gf4_poly_get_coefficient();
void test_gf4_poly_set_coefficient();
void test_gf4_poly_add();
void test_gf4_poly_add_inplace();
void test_gf4_poly_add_ax_to_deg_inplace();
void test_gf4_poly_mul();
void test_gf4_poly_div_x_to_deg();
void test_gf4_poly_div_x_to_deg_inplace();
void test_gf4_poly_div_rem();
void test_gf4_poly_invert_slow();
void test_gf4_poly_is_zero();
void test_gf4_poly_equal();
void test_gf4_poly_cyclic_shift_right_inplace();
void test_gf4_poly_get_degree();
void test_gf4_poly_adjust_degree();
void test_gf4_poly_clone();
void test_gf4_poly_copy();

// gf4_matrix
void test_gf4_square_matrix_init_cyclic_matrix();
void test_gf4_matrix_gaussian_elimination_inplace();
void test_gf4_matrix_solve_homogenous_linear_system();

// contexts
void test_contexts_init();
void test_contexts_save_load();

// enc
void test_enc_encode();
void test_enc_encrypt();

// dec
void test_dec_calculate_syndrome();


// test runner
void run_unit_tests();

#endif