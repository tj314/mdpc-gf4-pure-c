#ifndef MDPC_GF4_DEC_H
#define MDPC_GF4_DEC_H

#include "gf4.h"
#include "gf4_poly.h"
#include "contexts.h"

/**
 * Calculate the syndrome of a vector.
 *
 * out_syndrome must be initialized in advance. It must have capacity of at least block_size. It must be a zero polynomial.
 *
 * @param out_syndrome pointer to a polynomial that will be used to store the calculated syndrome
 * @param in_vector pointer to a polynomial representing a vector
 * @param ctx a valid decoding context
 */
void dec_calculate_syndrome(gf4_poly_t * out_syndrome, gf4_poly_t * in_vector, decoding_context_t * ctx);

long dec_calculate_threshold_0(long syndrome_weight);

long dec_calculate_threshold_1(long syndrome_weight);

long dec_calculate_threshold_2(long syndrome_weight);

long dec_calculate_threshold_3(long syndrome_weight);

long dec_calculate_threshold_4(long syndrome_weight);

long dec_calculate_threshold_5(long syndrome_weight);

long dec_calculate_new_sigma(gf4_poly_t * h_block, gf4_poly_t * syndrome, long syndrome_weight, gf4_t a, size_t actual_j, decoding_context_t * ctx);

void dec_flip_symbol(gf4_poly_t * h_block, gf4_poly_t * syndrome, gf4_poly_t * maybe_decoded, gf4_t a, size_t actual_j,size_t j, decoding_context_t * ctx);


/**
 * Perform basic symbol-flipping decoding.
 *
 * maybe_decoded and in_vector must be initialized in advance and must have capacity at least 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * @param maybe_decoded pointer to a polynomial that will be used to store the decoded message
 * @param in_vector pointer to a polynomial representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx);

/**
 * Perform basic symbol-flipping decoding, version 2.
 *
 * maybe_decoded and in_vector must be initialized in advance and must have capacity at least 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * while symbol flipping thies to maximize w(s) - w(s-a*h_j), this version tries to minimize w(s-a*h_j).
 * This version of SF decoder is DEPRECATED and will be removed.
 *
 * @param maybe_decoded pointer to a polynomial that will be used to store the decoded message
 * @param in_vector pointer to a polynomial representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping_2(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx);

/**
 * Perform symbol-flipping decoding with delta param.
 *
 * maybe_decoded and in_vector must be initialized in advance and must have capacity = 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * Symbol flipping searches maximum positive value of maximize w(s) - w(s-a*h_j) = sigma_max.
 * This version finds sigma_max and then flips all positions where sigma_j >= sigma_max-DELTA.
 * The required value of DELTA must be set in ctx before calling thus decoder!
 * It is recommended to use DELTA values s.t. 0 <= DELTA <= 3.
 * For instance, setting the value to 3: ctx->delta_setting = 3
 *
 * @param maybe_decoded pointer to a polynomial that will be used to store the decoded message
 * @param in_vector pointer to a polynomial representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping_delta(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx);

/**
 * Perform symbol-flipping decoding with delta param.
 *
 * maybe_decoded and in_vector must be initialized in advance and must have capacity = 2*block_size.
 * ctx must a valid decoding_context_t.
 *
 * This decoder calculates threshold T and then flips all positions where sigma_j > T.
 * The function to calculate the threshold value must be set in the decoding context prior to calling this function!
 * For instance, ctx->threshold = &dec_calculate_threshold_1;
 *
 * @param maybe_decoded pointer to a polynomial that will be used to store the decoded message
 * @param in_vector pointer to a polynomial representing the encoded message
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decoding, false otherwise
 */
bool dec_decode_symbol_flipping_threshold(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx);

/**
 * Decrypt an encrypted message using specified decoder.
 *
 * out_decrypted and in_encrypted must be initialized in advance and must have capacity at least 2*ctx->block size!
 * After a successful decryption, the decrypted messige is store as the first ctx->block_size symbols in out_encrypted.
 * The remaining ctx->block_size symbols are set to 0. out_decrypted->degree MAY NOT be set!
 * In the case of unsuccessful decryption (decoding failure was encountered), all symbols in out_decrypted are set to 0.
 *
 * @param out_decrypted pointer a polynomial that will be used to store the decrypted message
 * @param in_encrypted pointer to a polynomial representing the ciphertext
 * @param decode pointer to the decoder function to be used
 * @param num_iterations number of decoding iterations
 * @param ctx a valid decoding context
 * @return true on successful decryption, false otherwise
 */
bool dec_decrypt(gf4_poly_t * out_decrypted, gf4_poly_t * in_encrypted, bool (*decode)(gf4_poly_t*, gf4_poly_t *, size_t, decoding_context_t *), size_t num_iterations, decoding_context_t * ctx);


bool dec_decode_symbol_flipping_bgf(gf4_poly_t * maybe_decoded, gf4_poly_t * in_vector, size_t num_iterations, decoding_context_t * ctx);
#endif //MDPC_GF4_DEC_H
