#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "src/contexts.h"
#include "src/utils.h"
#include "src/random.h"
#include "src/dec.h"

const size_t block_size = 2293; // 2339;
const size_t num_errors = 84; // 96;
size_t M = 500;
size_t dist_start;
size_t dist_stop;
bool same;

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
        contexts_init(&ec, &dc, block_size, 37);
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

void print_usage() {
    fprintf(stderr,
            "./mdpc-gf4 gen\n"
            "./mdpc-gf4 same <dist-start> <dist-stop> <M>\n"
            "./mdpc-gf4 diff <dist-start> <dist-stop> <M>\n"
            "\n"
            "dist-start is inclusive\n"
            "dist-stop is exclusive\n"
            "M is number of messages per distance\n"
            "example: ./mdpc-gf4 same 1 3 1000 --> use 1000 error patterns targeting same symbols for distances 1 and 2\n");
}

void syndrome_weights_sums_per_distance(const char * fname, void (*errvec_gen)(gf4_poly_t*, size_t, size_t, size_t), decoding_context_t * dc) {
    assert(NULL != fname);
    assert(dist_start < dist_stop);
    gf4_poly_t errvec = gf4_poly_init_zero(2*block_size);
    gf4_poly_t synvec = gf4_poly_init_zero(block_size);

    size_t * syn_weight_sums = malloc((dist_stop - dist_start) * sizeof(size_t *));
    assert(NULL != syn_weight_sums);

    size_t end = dist_stop - dist_start;
    for (size_t dist = dist_start; dist < dist_stop; ++dist) {
        size_t idx = dist - dist_start;
        fprintf(stderr, "Progress: %04zu / %04zu   ", idx + 1, end);
        for (size_t i = 0; i < M; ++i) {
            errvec_gen(&errvec, block_size, num_errors, dist);
            dec_calculate_syndrome(&synvec, &errvec, dc);

            syn_weight_sums[idx] += utils_hamming_weight(&synvec);
            gf4_poly_zero_out(&synvec);
            gf4_poly_zero_out(&errvec);
        }
        fprintf(stderr, "sum = %zu\n", syn_weight_sums[idx]);
    }
    FILE * out = fopen(fname, "w");
    for (size_t dist = dist_start; dist < dist_stop; ++dist) {
        fprintf(out, "%zu:%zu\n", dist, syn_weight_sums[dist - dist_start]);
    }
    fclose(out);
    gf4_poly_deinit(&errvec);
    gf4_poly_deinit(&synvec);
    free(syn_weight_sums);
}

int main(int nargs, char ** argv) {
    void (*f)(gf4_poly_t*, size_t, size_t, size_t) = NULL;
    if (2 == nargs && 0 == strcmp(argv[1], "gen")) {
        gen_keys("keys.txt", "mults.txt");
        return 0;
    } else if (5 == nargs && 0 == strcmp(argv[1], "same")) {
        f = &random_weighted_gf4_poly_pairs_of_ones;
    } else if (5 == nargs && 0 == strcmp(argv[1], "one-alpha")) {
        f = &random_weighted_gf4_poly_pairs_of_one_alpha;
    } else if (5 == nargs && 0 == strcmp(argv[1], "alpha-one")) {
        f = &random_weighted_gf4_poly_pairs_of_alpha_one;
    } else {
        print_usage();
        return 0;
    }

    dist_start = atoll(argv[2]);
    dist_stop = atoll(argv[3]);
    M = atoll(argv[4]);

    fprintf(stderr, "used settings: same=%hhu, dist-start=%zu, dist-stop=%zu, M=%zu\n", (uint8_t)same, dist_start, dist_stop, M);
    if (dist_stop <= dist_start) {
        fprintf(stderr, "ERROR: expected dist-start < dist-stop! Got dist-start=%zu and dist-stop=%zu!\n", dist_start, dist_stop);
        return -1;
    }

    encoding_context_t ec;
    decoding_context_t dc;
    contexts_load("keys.txt", block_size, &ec, &dc);

    char filename[100] = {0};
    sprintf(filename, "syndrome-weights-sums-%s-from-%zu-to-%zu.txt", argv[1], dist_start, dist_stop);

    syndrome_weights_sums_per_distance(filename, f, &dc);

    contexts_deinit(&ec, &dc);
    return 0;
}