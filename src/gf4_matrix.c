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
    matrix.rows = malloc(N * sizeof(gf4_t *));
    if (NULL == matrix.rows) {
        fprintf(stderr, "gf4_square_matrix_make_cyclic_matrix: Memory allocation failed!\n");
        exit(-1);
    }
    for (size_t row = N; row > 0; --row) {
        matrix.rows[N - row] = malloc(N * sizeof(gf4_t));
        if (NULL == matrix.rows[N - row]) {
            fprintf(stderr, "gf4_square_matrix_make_cyclic_matrix: Memory allocation failed!\n");
            exit(-1);
        }
        for (size_t col = 0; col < N; ++col) {
            matrix.rows[N - row][col] = first_row->coefficients[(row + col) % N];
        }
    }
    matrix.N = N;
    return matrix;
}

void gf4_square_matrix_deinit(gf4_square_matrix_t * matrix) {
    assert(NULL != matrix);
    for (size_t i = 0; i < matrix->N; ++i) {
        free(matrix->rows[i]);
        matrix->rows[i] = NULL;
    }
    free(matrix->rows);
    matrix->rows = NULL;
}