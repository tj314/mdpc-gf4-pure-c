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

gf4_matrix_t gf4_matrix_init_cyclic_matrix(gf4_poly_t * first_row, size_t N) {
    assert(NULL != first_row);
    assert(first_row->capacity >= N);
    assert(0 != N);
    gf4_matrix_t matrix;
    matrix.rows = malloc(N * sizeof(gf4_t *));
    if (NULL == matrix.rows) {
        fprintf(stderr, "gf4_matrix_init_cyclic_matrix: Memory allocation failed!\n");
        exit(-1);
    }
    for (size_t row = N; row > 0; --row) {
        matrix.rows[N - row] = malloc(N * sizeof(gf4_t));
        if (NULL == matrix.rows[N - row]) {
            fprintf(stderr, "gf4_matrix_init_cyclic_matrix: Memory allocation failed!\n");
            exit(-1);
        }
        for (size_t col = 0; col < N; ++col) {
            matrix.rows[N - row][col] = first_row->coefficients[(row + col) % N];
        }
    }
    matrix.num_rows = N;
    matrix.num_cols = N;
    return matrix;
}

bool gf4_matrix_find_pivot_index(gf4_matrix_t * matrix, size_t row_index, size_t * out_pivot_index) {
    // find the index of pivot in the given row
    for (size_t col = 0; col < matrix->num_cols; ++col) {
        if (0 != matrix->rows[row_index][col]) {
            *out_pivot_index = col;
            return true;
        }
    }
    *out_pivot_index = matrix->num_cols;
    return false;
}

bool gf4_matrix_is_row_zero(gf4_matrix_t * matrix, size_t row_index) {
    bool is_zero = true;
    for (size_t i = 0; i < matrix->num_rows; ++i) {
        is_zero &= (0 == matrix->rows[row_index][i]);
    }
    return is_zero;
}

bool gf4_matrix_is_upper_echelon(gf4_matrix_t * matrix) {
    assert(NULL != matrix);
    size_t pivot_index, last_pivot_index = 0;
    bool last_pivot_exists = true;
    for (size_t row = 0; row < matrix->num_rows; ++row) {
        bool pivot_exists = gf4_matrix_find_pivot_index(matrix, row, &pivot_index);
        if (last_pivot_exists) {
            last_pivot_exists = pivot_exists;
            if (last_pivot_index > pivot_index) {
                return false;
            }
            for (size_t current_row = 0; current_row < matrix->num_rows; ++current_row) {
                if (row == current_row) {
                    continue;
                }
                if (matrix->rows[current_row][pivot_index] != 0) {
                    return false;
                }
            }
            last_pivot_index = pivot_index;
        } else if (pivot_exists) {
            // last row was all zeros, but this row is not
            return false;
        }
    }
    return true;
}

/**
 * @brief find the index of minimum and invalidate the minimum
 * @param array
 * @param length
 * @param max_value
 * @return
 */
size_t gf4_matrix_helper_find_min(size_t * array, size_t length, size_t max_value) {
    size_t min = max_value;
    size_t min_index = length;
    for (size_t i = 0; i < length; ++i) {
        if (array[i] < min) {
            min_index = i;
            min = array[i];
        }
    }
    if (min_index < length) {
        array[min_index] = max_value;
    }
    return min_index;
}

void gf4_matrix_gaussian_elimination_inplace(gf4_matrix_t * matrix) {
    size_t * pivot_indices = calloc(matrix->num_rows, sizeof(size_t));
    assert(NULL != pivot_indices);
    // apply gaussian elimination from top to bottom
    for (size_t row = 0; row < matrix->num_rows; ++row) {
        size_t pivot_index;
        bool pivot_found = gf4_matrix_find_pivot_index(matrix, row, &pivot_index);
        pivot_indices[row] = pivot_index;
        if (!pivot_found) {
            // zero row
            continue;
        }
        for (size_t other_row = 0; other_row < matrix->num_rows; ++other_row) {
            if (row == other_row || 0 == matrix->rows[other_row][pivot_index]) {
                continue;
            }
            gf4_t ratio = gf4_div(matrix->rows[row][pivot_index], matrix->rows[other_row][pivot_index]);
            for (size_t col = 0; col < matrix->num_cols; ++ col) {
                matrix->rows[other_row][col] = matrix->rows[row][col] ^ gf4_mul(matrix->rows[other_row][col], ratio);
            }
        }
    }

    // swap rows to obtain echelon form
    size_t min_index = gf4_matrix_helper_find_min(pivot_indices, matrix->num_rows, matrix->num_cols);
    size_t row_pos = 0;
    while (min_index < matrix->num_rows) {
        // swap rows in matrix
        gf4_t * tmp_row = matrix->rows[row_pos];
        matrix->rows[row_pos] = matrix->rows[min_index];
        matrix->rows[min_index] = tmp_row;
        // swap pivot indices
        size_t tmp_index = pivot_indices[row_pos];
        pivot_indices[row_pos] = pivot_indices[min_index];
        pivot_indices[min_index] = tmp_index;
        // update row_pos and min_index
        min_index = gf4_matrix_helper_find_min(pivot_indices, matrix->num_rows, matrix->num_cols);
        ++row_pos;
    }

    // cleanup
    free(pivot_indices);
}

size_t gf4_matrix_rank(gf4_matrix_t * matrix) {
    size_t num_zero = 0;
    for (size_t row = 0; row < matrix->num_rows; ++row) {
        num_zero += (size_t)gf4_matrix_is_row_zero(matrix, row);
    }
    assert(matrix->num_rows >= num_zero);
    return matrix->num_rows - num_zero;
}

gf4_matrix_t gf4_matrix_solve_homogenous_linear_system(gf4_matrix_t * equations) {
    assert(NULL != equations);
    assert(gf4_matrix_is_upper_echelon(equations));
    size_t equations_rank = gf4_matrix_rank(equations);
    size_t num_params = equations->num_cols - equations_rank;

    // allocate matrix for the solutions
    gf4_matrix_t out_solutions;
    out_solutions.num_cols = equations->num_cols;
    out_solutions.num_rows = utils_binary_pow(4, num_params);
    out_solutions.rows = malloc(out_solutions.num_rows * sizeof(gf4_t *));
    assert(NULL != out_solutions.rows);
    for (size_t i = 0; i < out_solutions.num_rows; ++i) {
        out_solutions.rows[i] = calloc(out_solutions.num_cols, sizeof(gf4_t));
        assert(NULL != out_solutions.rows[i]);
    }

    // solve
    if (0 == equations_rank) {
        fprintf(stderr, "gf4_matrix_solve_homogenous_linear_system: The solution to the system of equations is all permutations of length %zu. This is not yet implemented!\n", equations->num_cols);
        exit(-1);
    } else {
        size_t current_row = equations_rank;
        size_t repetitions = 1;
        size_t last_set_unknown_index = equations->num_cols - 1;
        do {
            --current_row;
            size_t pivot_index;
            bool pivot_exists = gf4_matrix_find_pivot_index(equations, current_row, &pivot_index);
            if (!pivot_exists) {
                fprintf(stderr, "gf4_matrix_solve_homogenous_linear_system: Unexpected error! Encountered a zero row where no such row was expected!\n");
                exit(-1);
            }
            for (size_t col = last_set_unknown_index; col > pivot_index; --col) {
                gf4_t val = 0;
                size_t index = 0;
                while (index < out_solutions.num_rows) {
                    for (size_t rep = 0; rep < repetitions; ++rep) {
                        out_solutions.rows[index][col] = val;
                        ++index;
                    }
                    ++val;
                    val %= 4;
                }
                repetitions *= 4;
            }
            for (size_t row = 0; row < out_solutions.num_rows; ++row) {
                gf4_t sum_right_of_pivot = 0;
                for (size_t col = pivot_index + 1; col < out_solutions.num_cols; ++col) {
                    sum_right_of_pivot ^= gf4_mul(equations->rows[row][col], out_solutions.rows[row][col]);
                }
                out_solutions.rows[row][pivot_index] = gf4_div(sum_right_of_pivot, equations->rows[row][pivot_index]);
            }
            last_set_unknown_index = (0 == pivot_index) ? pivot_index - 1 : out_solutions.num_cols;
        } while (0 != current_row);
    }


    return out_solutions;
}

void gf4_matrix_deinit(gf4_matrix_t * matrix) {
    assert(NULL != matrix);
    for (size_t i = 0; i < matrix->num_rows; ++i) {
        free(matrix->rows[i]);
        matrix->rows[i] = NULL;
    }
    free(matrix->rows);
    matrix->rows = NULL;
}

void gf4_matrix_pretty_print(gf4_matrix_t * matrix) {
    for (size_t row = 0; row < matrix->num_rows; ++row) {
        for (size_t col = 0; col < matrix->num_cols; ++col) {
            printf("%5s ", gf4_to_str(matrix->rows[row][col]));
        }
        printf("\n");
    }
}