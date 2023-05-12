#ifndef MDPC_GF4_ENC_H
#define MDPC_GF4_ENC_H

#include "gf4.h"
#include "gf4_poly.h"
#include "contexts.h"

/**
 * @brief Encode a message.
 *
 * out_encoded and in_message must be initialized in advance with capacity at least 2*ctx->block_size and ctx->block_size respectively.
 * out_encoded must be a zero polynomial.
 *
 * @param out_encoded vector to store the result in
 * @param in_message vector containing a message to encode
 * @param ctx a valid encoding context
 */
void enc_encode(gf4_poly_t * out_encoded, gf4_poly_t * in_message, encoding_context_t * ctx);

/**
 * @brief Encrypt a message.
 *
 * out_encrypted and in_message must be initialized in advance with capacity at least 2*ctx->block_size and ctx->block_size respectively.
 * out_encrypted must be a zero polynomial.
 * in_message is first encoded. A random error vector with hamming weight equal to num_errors is then generated.
 * The encrypted message is then equal to encoded+error_vector.
 *
 * @param out_encrypted vector to store the result in
 * @param in_message vector containing a message to encrypt
 * @param num_errors hamming weight of the error vector to be used
 * @param ctx a valid encoding context
 */
void enc_encrypt(gf4_poly_t * out_encrypted, gf4_poly_t * in_message, size_t num_errors, encoding_context_t * ctx);

#endif //MDPC_GF4_ENC_H
