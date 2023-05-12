/**
 *  @file   contexts.h
 *  @brief  Encoding and decoding contexts.
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


#ifndef MDPC_GF4_CONTEXTS_H
#define MDPC_GF4_CONTEXTS_H

#include <stdlib.h>
#include <inttypes.h>
#include "gf4.h"
#include "gf4_poly.h"
#include "utils.h"
#include "random.h"

/**
 * Structure that holds the public key.
 *
 * This structure holds the second part of matrix G.
 * In general, G = (I | -(H0*H1^-1)^T) (where H=(H0|H1))
 * Here, second_block_G represents the second block of matrix G,
 * but it is not transposed: second_block_G = H0*H1^-1.
 * Transposition is not necessary, the encoding function was written
 * to accommodate for this.
 *
 */
typedef struct {
    gf4_poly_t second_block_G; ///< polynomial representing the first row of the second block of matrix G, not transposed
    size_t block_size; ///< size of the circulant block
#ifdef WRITE_WEIGHTS
        size_t index; ///< index to distinguish various runs of experiments
#endif
} encoding_context_t;

/**
 * @brief Structure that holds the private key.
 *
 * This structure holds the matrix H.
 * In general: H = (H0 | H1 | ... | HN)
 * In this implementation, matrix H consists of two circulant blocks,
 * as is customary in related cryptosystems such as BIKE, i.e. H = (H0 | H1).
 *
 */
typedef struct {
    gf4_poly_t h0; ///< polynomial representing the first row of the first block of matrix H
    gf4_poly_t h1; ///< polynomial representing the first row of the second block of matrix H
    size_t block_size; ///< size of the circulant block
    long delta_setting; ///< setting for the parameter delta used by some decoders
    long (*threshold)(long); ///< function to calculate the threshold based on syndrome weight used by some decoders
#ifdef WRITE_WEIGHTS
    size_t index; ///< index to distinguish various runs of experiments
#endif
    size_t elapsed_iterations; ///< number of elapsed iterations during decoding, may be set by some decoders
} decoding_context_t;

/**
 * @brief Generate contexts for encoding and decoding.
 *
 * generate polynomials h0 and h1 st. hamming weight of h0 (and also h1) is equal to block_weight.
 * generated polynomial h1 will be invertible mod (x^block_size + 1). h0 and h1 are used for decoding.
 * polynomial for encoding is calculated as follows: (h0 * h1^-1) mod (x^block_size + 1).
 * This function will allocate memory for the polynomials. Do not initialize polynomials in out_enc_ctx and
 * out_dec_ctx yourself!
 * Members of encoding_context_t and decoding_context_t that are not relevant to all decoders may not be initialized!
 *
 * @see contexts_deinit
 *
 * @param out_enc_ctx memory location to store encoding polynomials and parameters to
 * @param out_dec_ctx memory location to store decoding polynomials and parameters to
 * @param block_size size of the circulant block of the matrices H, G
 * @param block_weight hamming weight of each row/columns of the circulant blocks of H and G
 */
void contexts_init(encoding_context_t * out_enc_ctx, decoding_context_t * out_dec_ctx, size_t block_size, size_t block_weight);

/**
 * @brief Deinit contexts.
 *
 * Deinitializes polynomials contained within the contexts.
 *
 * @param enc_ctx memory location of the encoding context
 * @param dec_ctx memory location of the decoding context
 */
void contexts_deinit(encoding_context_t * enc_ctx, decoding_context_t * dec_ctx);

/**
 * @brief Save matrices G and H to a file.
 *
 * @param filename savefile path
 * @param enc_ctx memory location of the encoding context
 * @param dec_ctx memory location of the decoding context
 */
void contexts_save(const char * filename, encoding_context_t * enc_ctx, decoding_context_t * dec_ctx);

/**
 * @brief Load matrices G and H to a file.
 *
 * Allocates all the necessary memory for enc_ctx and dec_ctx. Do not allocate them yourself!
 * Do not call contexts_init with the same enc_ctx and dec_ctx!
 *
 * @see contexts_deinit
 *
 * @param filename savefile path
 * @param enc_ctx memory location of the encoding context
 * @param dec_ctx memory location of the decoding context
 */
void contexts_load(const char * filename, size_t block_size, encoding_context_t * enc_ctx, decoding_context_t * dec_ctx);
#endif // MDPC_GF4_CONTEXTS_H
