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

#ifndef GF4_GF4_MATRIX_H
#define GF4_GF4_MATRIX_H

#include "gf4_poly.h"
#include "utils.h"

/**
 * @brief This structure represents a square matrix of size NxN.
 */
typedef struct {
    gf4_t ** rows; ///< rows of the matrix
    size_t num_rows; ///< number of rows
    size_t num_cols; ///< number of columns
} gf4_matrix_t;

/**
 * @brief Allocate a cyclic matrix and fill it it based on the given first_row.
 *
 * i-th row of the cyclic matrix is obtained by cyclically shifting its (i-1)-th row by 1 place to the right.
 * This function allocates memory!
 * Initialized matrix must be cleaned up using gf4_matrix_deinit function if no longer needed!
 *
 * @see gf4_matrix_deinit
 *
 * @param first_row the first row of the matrix
 * @param N size of the matrix
 * @return initialized cyclic matrix constructed from the provided first row
 */
gf4_matrix_t gf4_matrix_init_cyclic_matrix(gf4_poly_t * first_row, size_t N);

/**
 * @brief Destroy a matrix.
 *
 * Free all the rows of the matrix and the rows array itself.
 *
 * @param matrix an initialized matrix
 */
void gf4_matrix_deinit(gf4_matrix_t * matrix);


/**
 * @brief Perform gaussian elimination on the matrix.
 *
 * The resulting matrix will be in row echelon form (all values above and bellow each row's pivot are 0).
 *
 * @param matrix an initialized matrix
 */
void gf4_matrix_gaussian_elimination_inplace(gf4_matrix_t * matrix);


/**
 * @brief Solve a system of linear homogenous equations.
 *
 * The matrix of equations is expected to be in row echelon form.
 * The solutions are stored in a new matrix with 4^(equations->num_columns - rank) rows and equations->num_columns columns.
 * Therefore this function allocates memory!
 * Initialized matrix must be cleaned up using gf4_matrix_deinit function if no longer needed!
 *
 * @see gf4_matrix_gaussian_elimination_inplace
 * @see gf4_matrix_deinit
 *
 * @param equations an initialized matrix representing the system of linear equations in row echelon form
 * @return an initialized matrix of solutions
 */
gf4_matrix_t gf4_matrix_solve_homogenous_linear_system(gf4_matrix_t * equations);

/**
 * @brief Pretty print a matrix.
 *
 * @param matrix an initialized matrix
 */
void gf4_matrix_pretty_print(gf4_matrix_t * matrix);
#endif //GF4_GF4_MATRIX_H
