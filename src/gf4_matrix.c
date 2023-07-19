/**
 *  @file   qc_matrix.h
 *  @brief  QC (quasi-cyclic) matrix over GF(4).
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

#include "gf4_matrix.h"

gf4_square_matrix_t gf4_square_matrix_make_cyclic_matrix(gf4_poly_t * first_row, size_t N) {
    assert(NULL != first_row);
    assert(first_row->capacity >= N);
    assert(0 != N);
    gf4_square_matrix_t matrix;
    matrix.N = N;
    matrix.rows = calloc(N, sizeof(gf4_poly_t));
    if (NULL == matrix.rows) {
        fprintf(stderr, "gf4_square_matrix_make_cyclic_matrix: Memory allocation failed!\n");
        exit(-1);
    }
    matrix.rows[0] = gf4_poly_clone(first_row);
    for (size_t i = 1; i < N; ++i) {
        matrix.rows[i] = gf4_poly_clone(&matrix.rows[i-1]);
        gf4_poly_cyclic_shift_right_inplace(&matrix.rows[i], N);
    }
    return matrix;
}

void gf4_square_matrix_deinit(gf4_square_matrix_t * matrix) {
    assert(NULL != matrix);
    for (size_t i = 0; i < matrix->N; ++i) {
        if (NULL != &matrix->rows[i]) {
            gf4_poly_deinit(&matrix->rows[i]);
        }
    }
    free(matrix->rows);
    matrix->rows = NULL;
}