#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "src/contexts.h"
#include "src/utils.h"
#include "src/random.h"
#include "src/dec.h"
#include "src/gf4_matrix.h"

#define IS_ALPHA_MULT_RIGHT(a, b) ((1 == a && 2 == b) || (2 == a && 3 == b) || (3 == a && 1 == b))
#define IS_ALPHA_MULT_LEFT(a, b) ((2 == a && 1 == b) || (3 == a && 2 == b) || (1 == a && 3 == b))
#define ABS_DIFF(a, b) ((a > b) ? a - (b) : b - (a))

const size_t block_size =  2339; // 2293;
const size_t num_errors = 84; // 96;
const size_t block_weight = 37;
size_t M = 1000;

// compare two size_t values; used for qsort
int compare_func(const void * a, const void * b) {
    size_t aa = *(size_t *)a;
    size_t bb = *(size_t *)b;
    if (aa > bb) {
        return 1;
    } else if (aa == bb) {
        return 0;
    } else {
        return -1;
    }
}

void gen_keys(const char * keys_fname, const char * mults_fname) {
    size_t ** mults_same, ** mults_diff;
    mults_same = malloc(4 * sizeof(size_t *));
    mults_diff = malloc(4 * sizeof(size_t *));
    assert(NULL != mults_same);
    assert(NULL != mults_diff);
    for (size_t i = 0; i < 4; ++i) {
        mults_same[i] = calloc(block_size, sizeof(size_t));
        mults_diff[i] = calloc(block_size, sizeof(size_t));
        assert(NULL != mults_same[i]);
        assert(NULL != mults_diff[i]);
    }
    bool end = false;
    while (true) {
        encoding_context_t ec;
        decoding_context_t dc;
        contexts_init(&ec, &dc, block_size, block_weight);
        utils_get_distance_multiplicities_h0(mults_same, mults_diff, &dc);
        for (size_t i = 0; i < 3; ++i) {  // 3 is correct
            for (size_t j = 0; j < 11; ++j) {
                if (0 == mults_same[i][j] || 0 == mults_diff[i][j]) {
                    goto REGEN;
                }
            }
        }
        end = true;
        contexts_save(keys_fname, &ec, &dc);
        FILE * mults_file = fopen(mults_fname, "w");
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 11; ++j) {
                fprintf(mults_file, "%zu ", mults_same[i][j]);
            }
            fprintf(mults_file, "\n");
        }
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 11; ++j) {
                fprintf(mults_file, "%zu ", mults_diff[i][j]);
            }
            fprintf(mults_file, "\n");
        }
        fclose(mults_file);
        REGEN:
        for (size_t i = 0; i < 4; ++i) {
            memset(mults_same[i], 0, sizeof(size_t) * block_size);
            memset(mults_diff[i], 0, sizeof(size_t) * block_size);
        }
        contexts_deinit(&ec, &dc);
        if (end) {
            break;
        }
    }
    for (size_t i = 0; i < 4; ++i) {
        free(mults_diff[i]);
        free(mults_same[i]);
    }
    free(mults_same);
    free(mults_diff);
}

void collect_syndrome_weights(const size_t id) {
    size_t half_block_size = (block_size / 2) + 1; // half of the block size rounded up

    size_t * weights_distances_same_e0 = calloc(half_block_size, sizeof(size_t));
    size_t * weights_distances_alpha_multiple_right_e0 = calloc(half_block_size, sizeof(size_t));
    size_t * weights_distances_alpha_multiple_left_e0 = calloc(half_block_size, sizeof(size_t));

    size_t * weights_distances_same_e1 = calloc(half_block_size, sizeof(size_t));
    size_t * weights_distances_alpha_multiple_right_e1 = calloc(half_block_size, sizeof(size_t));
    size_t * weights_distances_alpha_multiple_left_e1 = calloc(half_block_size, sizeof(size_t));

    if (NULL == weights_distances_same_e0
        || NULL == weights_distances_alpha_multiple_right_e0
        || NULL == weights_distances_alpha_multiple_left_e0
        || NULL == weights_distances_same_e1
        || NULL == weights_distances_alpha_multiple_right_e1
        || NULL == weights_distances_alpha_multiple_left_e1) {
        fprintf(stderr, "%s (%d): Allocation error!\n", __func__, __LINE__);
        exit(-1);
    }

    gf4_array_t err = gf4_array_init(2*block_size, true);
    gf4_array_t syndrome = gf4_array_init(block_size, true);

    encoding_context_t ec;
    decoding_context_t dc;
    contexts_load("keys.txt", block_size, &ec, &dc);

    for (size_t run = 0; run < M; ++run) {
        fprintf(stderr, "Progress: %zu/%zu\n", (run+1), M);
        random_weighted_gf4_array(&err, 2*block_size, num_errors);
        dec_calculate_syndrome(&syndrome, &err, &dc);
        size_t syndrome_weight = gf4_array_hamming_weight(&syndrome);

        for (size_t i = 0; i < block_size; ++i) {
            for (size_t j = i + 1; j < block_size; ++j) {
                size_t distance = (j - i) < half_block_size ? j - i : block_size - (j - i);

                // first half of err (e0)
                gf4_t a_e0 = err.array[i], b_e0 = err.array[j];
                if (0 != a_e0 && a_e0 == b_e0) {
                    weights_distances_same_e0[distance] += syndrome_weight;
                } else if (IS_ALPHA_MULT_RIGHT(a_e0, b_e0)) {
                    weights_distances_alpha_multiple_right_e0[distance] += syndrome_weight;
                } else if (IS_ALPHA_MULT_LEFT(a_e0, b_e0)) {
                    weights_distances_alpha_multiple_left_e0[distance] += syndrome_weight;
                }

                // second half of err (e1)
                gf4_t a_e1 = err.array[block_size + i], b_e1 = err.array[block_size + j];
                if (0 != a_e1 && a_e1 == b_e1) {
                    weights_distances_same_e1[distance] += syndrome_weight;
                } else if (IS_ALPHA_MULT_RIGHT(a_e1, b_e1)) {
                    weights_distances_alpha_multiple_right_e1[distance] += syndrome_weight;
                } else if (IS_ALPHA_MULT_LEFT(a_e1, b_e1)) {
                    weights_distances_alpha_multiple_left_e1[distance] += syndrome_weight;
                }
            }
        }

        gf4_array_zero_out(&syndrome);
        gf4_array_zero_out(&err);
    }

    char buffer[100] = {0};

    // number of trials
    sprintf(buffer, "num_trials_%zu.txt", id);
    FILE * file = fopen(buffer, "w+");
    if (NULL == file) {
        fprintf(stderr, "%s (%d): Fopen error!\n", __func__, __LINE__);
        exit(-1);
    }
    fprintf(file, "%zu\n", M);
    fclose(file);

    // same symbols
    memset(buffer, 0, 100 * sizeof(char));
    sprintf(buffer, "weights_same_e0_%zu.txt", id);
    file = fopen(buffer, "w+");
    if (NULL == file) {
        fprintf(stderr, "%s (%d): Fopen error!\n", __func__, __LINE__);
        exit(-1);
    }
    for (size_t dist = 1; dist < half_block_size; ++dist) {
        fprintf(file, "%zu %zu\n", dist, weights_distances_same_e0[dist]);
    }
    fclose(file);

    memset(buffer, 0, 100 * sizeof(char));
    sprintf(buffer, "weights_same_e1_%zu.txt", id);
    file = fopen(buffer, "w+");
    if (NULL == file) {
        fprintf(stderr, "%s (%d): Fopen error!\n", __func__, __LINE__);
        exit(-1);
    }
    for (size_t dist = 1; dist < half_block_size; ++dist) {
        fprintf(file, "%zu %zu\n", dist, weights_distances_same_e1[dist]);
    }
    fclose(file);

    // alpha mult right
    memset(buffer, 0, 100 * sizeof(char));
    sprintf(buffer, "weights_alpha_mult_right_e0_%zu.txt", id);
    file = fopen(buffer, "w+");
    if (NULL == file) {
        fprintf(stderr, "%s (%d): Fopen error!\n", __func__, __LINE__);
        exit(-1);
    }
    for (size_t dist = 1; dist < half_block_size; ++dist) {
        fprintf(file, "%zu %zu\n", dist, weights_distances_alpha_multiple_right_e0[dist]);
    }
    fclose(file);

    memset(buffer, 0, 100 * sizeof(char));
    sprintf(buffer, "weights_alpha_mult_right_e1_%zu.txt", id);
    file = fopen(buffer, "w+");
    if (NULL == file) {
        fprintf(stderr, "%s (%d): Fopen error!\n", __func__, __LINE__);
        exit(-1);
    }
    for (size_t dist = 1; dist < half_block_size; ++dist) {
        fprintf(file, "%zu %zu\n", dist, weights_distances_alpha_multiple_right_e1[dist]);
    }
    fclose(file);

    // alpha multiple left
    memset(buffer, 0, 100 * sizeof(char));
    sprintf(buffer, "weights_alpha_mult_left_e0_%zu.txt", id);
    file = fopen(buffer, "w+");
    if (NULL == file) {
        fprintf(stderr, "%s (%d): Fopen error!\n", __func__, __LINE__);
        exit(-1);
    }
    for (size_t dist = 1; dist < half_block_size; ++dist) {
        fprintf(file, "%zu %zu\n", dist, weights_distances_alpha_multiple_left_e0[dist]);
    }
    fclose(file);

    memset(buffer, 0, 100 * sizeof(char));
    sprintf(buffer, "weights_alpha_mult_left_e1_%zu.txt", id);
    file = fopen(buffer, "w+");
    if (NULL == file) {
        fprintf(stderr, "%s (%d): Fopen error!\n", __func__, __LINE__);
        exit(-1);
    }
    for (size_t dist = 1; dist < half_block_size; ++dist) {
        fprintf(file, "%zu %zu\n", dist, weights_distances_alpha_multiple_left_e1[dist]);
    }
    fclose(file);

    contexts_deinit(&ec, &dc);
    gf4_array_deinit(&syndrome);
    gf4_array_deinit(&err);
    free(weights_distances_alpha_multiple_left_e1);
    free(weights_distances_alpha_multiple_right_e1);
    free(weights_distances_same_e1);
    free(weights_distances_alpha_multiple_left_e0);
    free(weights_distances_alpha_multiple_right_e0);
    free(weights_distances_same_e0);
}

// D0[i] == true => i \in D0. equiv. for D1
void reconstruct_private_key(bool * D0, bool * D1, size_t s0, size_t s1) {
    size_t half_block_size = (block_size / 2) + 1; // half of the block size rounded up
    bool * Z0 = calloc(block_size, sizeof(bool));
    bool * not_Z1 = calloc(block_size, sizeof(bool)); // this is Z1'
    if (NULL == not_Z1 || NULL == Z0) {
        fprintf(stderr, "%s: Allocation error!\n", __func__);
        exit(-1);
    }

    encoding_context_t ec;
    decoding_context_t dc;
    contexts_load("keys.txt", block_size, &ec, &dc);

    // construct Z0 and Z1' (called not_Z1 in the code)
    // we will use this to delete some rows and columns of matrix B later
    for (size_t i = 0; i < block_size; ++i) {
        // Z0 portion
        size_t distance1_Z0 = (i - 1) < half_block_size ? i - 1 : block_size - (i - 1);
        size_t abs_difference_Z0 = ABS_DIFF(s0 + 1, i);
        size_t distance2_Z0 = abs_difference_Z0 < half_block_size ? abs_difference_Z0 : block_size - abs_difference_Z0;
        if (D0[distance1_Z0] || D0[distance2_Z0]) {
            Z0[i] = true;
        }

        // not_Z1 portion
        size_t distance1_not_Z1 = distance1_Z0;
        size_t abs_difference_Z1 = ABS_DIFF(s1 + 1, i);
        size_t distance2_not_Z1 = abs_difference_Z1 < half_block_size ? abs_difference_Z1 : block_size - abs_difference_Z1;
        if (!D1[distance1_not_Z1] && !D1[distance2_not_Z1]) {
            not_Z1[i] = true;
        }
    }

    // matrix B is the transpose of the second block of G
    // however, transposition is not necessary here as second_block_G was never transposed in the first place
    // see encoding_context_t
    gf4_matrix_t B = gf4_matrix_init_cyclic_matrix(&ec.second_block_G, block_size);

    // keep only rows whose indices are in Z1'
    size_t idx = block_size;
    do {
        --idx;
        if (!not_Z1[idx]) {
            gf4_matrix_remove_row_inplace(&B, idx);
        }
    } while (0 != idx);

    for (size_t p = 0; p < block_size; ++p) {
        gf4_matrix_t B_prime = gf4_matrix_clone(&B);

        // keep only cols whose indices are in Z0
        idx = block_size;
        do {
            --idx;
            if (!Z0[idx]) {
                gf4_matrix_remove_col_inplace(&B_prime, idx);
            }
        } while (0 != idx);

        // solve
        gf4_matrix_gaussian_elimination_inplace(&B_prime);
        gf4_matrix_t kernel = gf4_matrix_solve_homogenous_linear_system(&B_prime);
        if (1 == kernel.num_rows) {
            gf4_array_t tmp_do_not_deinit;
            tmp_do_not_deinit.capacity = kernel.num_cols;
            tmp_do_not_deinit.array = kernel.rows[0];
            size_t num_nonzero = gf4_array_hamming_weight(&tmp_do_not_deinit);
            if (num_nonzero <= block_weight) {
                // this is the correct key

                // let's reconstruct h1
                gf4_array_t h1 = gf4_array_init(block_size, true);
                size_t col = 0;
                for (size_t i = 0; i < block_size; ++i) {
                    if (not_Z1[i]) {
                        h1.array[i] = kernel.rows[0][col];
                        ++col;
                    }
                }

                // and write the calculated h1 to a file
                char fname[100] = {0};
                sprintf(fname, "reconstructed_h1_p_%zu.txt", p);
                FILE * output = fopen(fname, "w+");
                if (NULL == output) {
                    fprintf(stderr, "%s: Output file couldn't be created!\n", __func__);
                    exit(-1);
                }
                fprintf(output, "%zu\n", num_nonzero);
                for (size_t i = 0; i < h1.capacity; ++i) {
                    if (0 != h1.array[i]) {
                        fprintf(output, "%zu %hhu\n", i, h1.array[i]);
                    }
                }
                fclose(output);
                // todo reconstruct h0
            }
        }

        // cyclically shift Z0 by one position to the right
        bool tmp = Z0[block_size - 1];
        for (size_t i = block_size - 1; i > 0; --i) {
            Z0[i] = Z0[i - 1];
        }
        Z0[0] = tmp;

        // cleanup
        gf4_matrix_deinit(&kernel);
        gf4_matrix_deinit(&B_prime);
    }

    gf4_matrix_deinit(&B);
    contexts_deinit(&ec, &dc);
    free(not_Z1);
    free(Z0);
}

void print_usage() {
    fprintf(stderr,
            "./mdpc-gf4 gen\n"
            "./mdpc-gf4 weights <M> <ID>\n"
            "\n"
            "M is number of messages per distance\n"
            "ID is used to identify resulting filenames\n"
            "example: ./mdpc-gf4 weights 1000 1 --> use 1000 error patterns and use 1 in the resulting files' names\n");
}

int main(int nargs, char ** argv) {
    if (2 == nargs && 0 == strcmp(argv[1], "gen")) {
        gen_keys("keys.txt", "mults.txt");
        return 0;
    } else if (4 == nargs && 0 == strcmp(argv[1], "weights")) {
        const size_t id = atoll(argv[3]);
        M = atoll(argv[2]);
        fprintf(stderr, "Used settings: M=%zu id=%zu\n", M, id);
        collect_syndrome_weights(id);
    } else {
        print_usage();
        return 0;
    }
    return 0;
}