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
#include "gf4_poly.h"
#include "random.h"
#include "utils.h"

// TESTS
// gf4
void test_gf4_is_in_range();
void test_gf4_to_str();
void test_gf4_add();
void test_gf4_mul();
void test_gf4_div();

// utils
void test_utils_hamming_weight();

// gf4_poly
void test_gf4_poly_cyclic_shift_right_inplace();

// contexts
void test_contexts_init();
void test_contexts_save_load();

// enc
void test_enc_encode();
void test_enc_encrypt();

// dec
void test_dec_calculate_syndrome();

#endif