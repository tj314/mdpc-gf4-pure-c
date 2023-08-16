/**
 *  @file   utils.h
 *  @brief  Additional utilities.
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

#ifndef MDPC_GF4_UTILS_H
#define MDPC_GF4_UTILS_H

#include <stdlib.h>
#include "contexts.h"

void utils_get_distance_multiplicities_h0(size_t ** multiplicities_same_symbols, size_t ** multiplicities_different_symbols, decoding_context_t * dc);

size_t utils_binary_pow(size_t x, size_t n);

#define UTILS_SUBTRACT_OR_ZERO(a, b) ((a >= b) ? (a - b) : 0)

#endif //MDPC_GF4_UTILS_H
