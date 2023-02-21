#ifndef MDPC_GF4_ENCODER_H
#define MDPC_GF4_ENCODER_H

#include "gf4.h"
#include "gf4_poly.h"
#include "contexts.h"

/**
 * Encode a message.
 *
 * out_encoded must be initialized in advance with sufficient capacity. It must be a zero polynomial.
 *
 * @param out_encoded polynomial to store the result in
 * @param in_message polynomial containing a message to encode
 * @param ctx a valid encoding context
 */
void encode(gf4_poly_t * out_encoded, gf4_poly_t * in_message, encoding_context_t * ctx);

#endif //MDPC_GF4_ENCODER_H
