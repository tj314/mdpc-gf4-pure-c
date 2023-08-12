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

#include "utils.h"

size_t utils_hamming_weight(gf4_vector_t *vector) {
    size_t weight = 0;
    for (size_t i = 0; i < vector->capacity; ++i) {
        weight += (0 != vector->array[i]);
    }
    return weight;
}

void utils_get_distance_multiplicities_h0(size_t ** multiplicities_same_symbols, size_t ** multiplicities_different_symbols, decoding_context_t * dc) {
    assert(NULL != multiplicities_same_symbols);
    assert(NULL != multiplicities_different_symbols);
    assert(NULL != dc);

    size_t * same_symbols_dists = calloc(dc->block_size, sizeof(size_t));
    assert(NULL != same_symbols_dists);
    size_t * diff_symbols_dists = calloc(dc->block_size, sizeof(size_t));
    assert(NULL != diff_symbols_dists);
    for (long i = 0; i < (long)dc->h0.capacity; ++i) {
        if (0 == dc->h0.array[i]) {
            continue;
        }
        for (long j = i + 1; j < (long)dc->h0.capacity; ++j) {
            if (0 == dc->h0.array[j]) {
                continue;
            }
            size_t dist_a = j - i;
            size_t dist_b = i + dc->block_size;
            dist_b -= j;
            size_t dist = (dist_a < dist_b) ? dist_a : dist_b;
            if (dc->h0.array[i] == dc->h0.array[j]) {
                same_symbols_dists[dist] += 1;
            } else {
                diff_symbols_dists[dist] += 1;
            }
        }
    }
    size_t last_index_same_symbols[4] = {0};
    size_t last_index_diff_symbols[4] = {0};
    for (size_t i = 1; i < dc->block_size; ++i) {  // i == 0 is not valid
        size_t current_same_symbols_dist = same_symbols_dists[i];
        size_t current_diff_symbols_dist = diff_symbols_dists[i];
        if (current_same_symbols_dist <= 3) {
            multiplicities_same_symbols[current_same_symbols_dist][last_index_same_symbols[current_same_symbols_dist]] = i;
            last_index_same_symbols[current_same_symbols_dist] += 1;
        }
        if (current_diff_symbols_dist <= 3) {
            multiplicities_different_symbols[current_diff_symbols_dist][last_index_diff_symbols[current_diff_symbols_dist]] = i;
            last_index_diff_symbols[current_diff_symbols_dist] += 1;
        }
    }
}

size_t utils_binary_pow(size_t x, size_t n) {
    if(n == 0) {
        return 1;
    } else {
        size_t t = utils_binary_pow(x, n / 2);
        t = t * t;
        return n % 2 ? x * t : t;
    }
}