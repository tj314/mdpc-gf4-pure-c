/**
 *  @file   dec.h
 *  @brief  Decoders and associated functions.
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

#ifndef MDPC_GF4_DEC_H
#define MDPC_GF4_DEC_H

#include "gf4.h"
#include "gf4_array.h"
#include "gf4_poly.h"
#include "contexts.h"
#include "utils.h"

/**
 * @brief Calculate the syndrome of an array.
 *
 * out_syndrome must be initialized in advance. It must have capacity of at least block_size. It must be a zero vector.
 *
 * @param out_syndrome pointer to an array that will be used to store the calculated syndrome
 * @param in_message pointer to an array
 * @param ctx a valid decoding context
 */
void dec_calculate_syndrome(gf4_array_t *out_syndrome, gf4_array_t *in_message, decoding_context_t * ctx);

/**
 * @brief Calculate adaptive threshold.
 *
 * The threshold here is calculated as max(floor(f0(syndrome_weight)), 0),
 * where f0(syndrome_weight) = 0.0248577875*syndrome_weight - 29.1143817.
 * f0 was found via linear regression.
 *
 * @param syndrome_weight the current hamming weight of the syndrome
 * @return calculated adaptive threshold.
 */
long dec_calculate_threshold_0(long syndrome_weight);

/**
 * @brief Calculate adaptive threshold.
 *
 * The threshold here is calculated as max(floor(f0(syndrome_weight) + 1), 0),
 * where f0(syndrome_weight) = 0.0248577875*syndrome_weight - 29.1143817.
 * f0 was found via linear regression.
 *
 * @param syndrome_weight the current hamming weight of the syndrome
 * @return calculated adaptive threshold.
 */
long dec_calculate_threshold_1(long syndrome_weight);

/**
 * @brief Calculate adaptive threshold.
 *
 * The threshold here is calculated as max(floor(f0(syndrome_weight) + 2), 0),
 * where f0(syndrome_weight) = 0.0248577875*syndrome_weight - 29.1143817.
 * f0 was found via linear regression.
 *
 * @param syndrome_weight the current hamming weight of the syndrome
 * @return calculated adaptive threshold.
 */
long dec_calculate_threshold_2(long syndrome_weight);

/**
 * @brief Calculate adaptive threshold.
 *
 * The threshold here is calculated as max(floor(f0(syndrome_weight) + 3), 0),
 * where f0(syndrome_weight) = 0.0248577875*syndrome_weight - 29.1143817.
 * f0 was found via linear regression.
 *
 * @param syndrome_weight the current hamming weight of the syndrome
 * @return calculated adaptive threshold.
 */
long dec_calculate_threshold_3(long syndrome_weight);

/**
 * @brief Calculate adaptive threshold.
 *
 * The threshold here is calculated as max(floor(f0(syndrome_weight) + 4), 0),
 * where f0(syndrome_weight) = 0.0248577875*syndrome_weight - 29.1143817.
 * f0 was found via linear regression.
 *
 * @param syndrome_weight the current hamming weight of the syndrome
 * @return calculated adaptive threshold.
 */
long dec_calculate_threshold_4(long syndrome_weight);

/**
 * @brief Calculate adaptive threshold.
 *
 * The threshold here is calculated as max(floor(f0(syndrome_weight) + 5), 0),
 * where f0(syndrome_weight) = 0.0248577875*syndrome_weight - 29.1143817.
 * f0 was found via linear regression.
 *
 * @param syndrome_weight the current hamming weight of the syndrome
 * @return calculated adaptive threshold.
 */
long dec_calculate_threshold_5(long syndrome_weight);

/**
 * @brief Calculate sigma_j.
 *
 * sigma_j = w(syndrome) - w(syndrome - a*H_j), where H_j is j-th column of matrix H.
 *
 * WARNING: This is an internal function for use within the decoders.
 * If you feel the need to use it in your client code, chances are, you are likely doing something very wrong.
 *
 * @param h_block pointer to the correct block in H
 * @param syndrome pointer to the allocated syndrome
 * @param syndrome_weight weight the syndrome
 * @param a gf4_t value, symbol to use
 * @param actual_j if h_block points to H0, actual_j = j, otherwise actual_j = j - cts->block_size
 * @param ctx a valid decoding context
 * @return calculated sigma_j
 */
long dec_calculate_new_sigma(gf4_poly_t * h_block, gf4_array_t *syndrome, long syndrome_weight, gf4_t a, size_t actual_j, decoding_context_t * ctx);

/**
 * @brief Flip the symbol in the output vector and update the syndrome.
 *
 * WARNING: This is an internal function for use within the decoders.
 * If you feel the need to use it in your client code, chances are, you are likely doing something very wrong.
 *
 * @param h_block pointer to the correct block in H
 * @param syndrome pointer to the allocated syndrome
 * @param maybe_decoded output vector to flip the symbol in
 * @param a gf4_t value, symbol to use
 * @param actual_j if h_block points to H0, actual_j = j, otherwise actual_j = j - cts->block_size
 * @param j position to flip
 * @param ctx a valid decoding context
 */
void dec_flip_symbol(gf4_poly_t * h_block, gf4_array_t *syndrome, gf4_array_t *maybe_decoded, gf4_t a, size_t actual_j, size_t j, decoding_context_t * ctx);


/**
 * @brief Perform basic symbol-flipping decoding.
 *
 * maybe_decoded and in_array must be initialized in advance and must have capacity at least 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * @param maybe_decoded pointer to an array that will be used to store the decoded message
 * @param in_array pointer to an array representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping(gf4_array_t *maybe_decoded, gf4_array_t *in_array, size_t num_iterations, decoding_context_t * ctx);

/**
 * @brief Perform symbol-flipping decoding with delta param.
 *
 * maybe_decoded and in_array must be initialized in advance and must have capacity = 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * Symbol flipping searches maximum positive value of maximize w(s) - w(s-a*h_j) = sigma_max.
 * This version finds sigma_max and then flips all positions where sigma_j >= sigma_max-DELTA.
 * The required value of DELTA must be set in ctx before calling this decoder!
 * For instance, setting the value to 3: ctx->delta_setting = 3
 *
 * @param maybe_decoded pointer to an array that will be used to store the decoded message
 * @param in_array pointer to an array representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping_delta(gf4_array_t *maybe_decoded, gf4_array_t *in_array, size_t num_iterations, decoding_context_t * ctx);

/**
 * @brief Perform symbol-flipping decoding with adaptive threshold.
 *
 * maybe_decoded and in_array must be initialized in advance and must have capacity = 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * This decoder calculates adaptive threshold T and then flips all positions where sigma_j > T.
 * The function to calculate the threshold value must be set in the decoding context prior to calling this function!
 * For instance, ctx->threshold = &dec_calculate_threshold_1;
 *
 * @param maybe_decoded pointer to an array that will be used to store the decoded message
 * @param in_array pointer to an array representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping_threshold(gf4_array_t *maybe_decoded, gf4_array_t *in_array, size_t num_iterations, decoding_context_t * ctx);

/**
 * @brief Perform symbol-flipping with black and gray symbols.
 *
 * maybe_decoded and in_array must be initialized in advance and must have capacity = 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * This decoder calculates adaptive threshold T and then flips all positions where sigma_j > T.
 * It saves the flipped symbols in the array called black.
 * It saves symbols that were almost flipped (sigma_j > T - DELTA) in the array called gray.
 * It then recalculates sigma_j for black symbols and flips some of them back.
 * It then recalculates sigma_j for gray symbols and flips some of them.
 * The function to calculate the threshold value must be set in the decoding context prior to calling this function!
 * For instance, ctx->threshold = &dec_calculate_threshold_1;
 *
 * WARNING: This function is under development and is not tested. Do not use it!
 *
 * @param maybe_decoded pointer to an array that will be used to store the decoded message
 * @param in_array pointer to an array representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping_bg(gf4_array_t *maybe_decoded, gf4_array_t *in_array, size_t num_iterations, decoding_context_t * ctx);

/**
 * @brief Decrypt an encrypted message using specified decoder.
 *
 * out_decrypted and in_encrypted must be initialized in advance and must have capacity at least 2*ctx->block size!
 * After a successful decryption, the decrypted messige is store as the first ctx->block_size symbols in out_encrypted.
 * The remaining ctx->block_size symbols are set to 0. out_decrypted->degree MAY NOT be set!
 * In the case of unsuccessful decryption (decoding failure was encountered), all symbols in out_decrypted are set to 0.
 *
 * @param out_decrypted pointer an array that will be used to store the decrypted message
 * @param in_encrypted pointer to an array representing the ciphertext
 * @param decode pointer to the decoder function to be used
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decryption, false otherwise
 */
bool dec_decrypt(gf4_array_t *out_decrypted, gf4_array_t *in_encrypted, bool (*decode)(gf4_array_t *, gf4_array_t *, size_t, decoding_context_t *), size_t num_iterations, decoding_context_t * ctx);
#endif //MDPC_GF4_DEC_H
