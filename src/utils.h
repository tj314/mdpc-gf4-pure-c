/**
 *  @file   utils.h
 *  @brief  Additional utilities.
 *  @author Tomáš Vavro
 *  @date   2023-05-12
 ***********************************************/

#ifndef MDPC_GF4_UTILS_H
#define MDPC_GF4_UTILS_H

#include <stdlib.h>
#include "gf4_poly.h"

/**
 * @brief Find Hamming weight of a vector.
 *
 * vector must be initialized beforehand.
 *
 * @param vector pointer to a vector
 * @return Hamming weight of vector
 */
size_t utils_hamming_weight(gf4_poly_t * vector);

#endif //MDPC_GF4_UTILS_H
