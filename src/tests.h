//#ifdef RUNTESTS
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

// contexts
void test_contexts_init();
void test_contexts_save_load();

// enc
void test_enc_encode();
void test_enc_encrypt();

// dec
void test_dec_calculate_syndrome();

#endif
//#endif