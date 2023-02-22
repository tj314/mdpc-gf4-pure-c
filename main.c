#include <stdio.h>
#include "src/gf4.h"
#include "src/gf4_poly.h"
#include "src/contexts.h"
#include "src/enc.h"
#include "src/dec.h"

int main(void) {
    /*
    gf4_poly_t poly1 = gf4_poly_init_zero(10);
    gf4_poly_t poly2 = gf4_poly_init_zero(10);
    gf4_poly_t modulus = gf4_poly_init_zero(10);
    gf4_poly_t maybe_inverse = gf4_poly_init_zero(10);

    gf4_poly_set_coefficient(&poly1, 2, 1);
    gf4_poly_set_coefficient(&poly1, 1, 1);
    gf4_poly_set_coefficient(&poly1, 0, 1);

    gf4_poly_set_coefficient(&poly2, 4, 2);
    gf4_poly_set_coefficient(&poly2, 1, 2);

    gf4_poly_set_coefficient(&modulus, 0, 1);
    gf4_poly_set_coefficient(&modulus, 8, 1);

    gf4_poly_t * poly = &poly1;
    bool inverse_exists = gf4_poly_invert_slow(&maybe_inverse, poly, &modulus);


    printf("poly:    ");
    gf4_poly_pretty_print(poly, "\n");
    printf("modulus: ");
    gf4_poly_pretty_print(&modulus, "\n");

    if (inverse_exists) {
        printf("inverse: ");
        gf4_poly_pretty_print(&maybe_inverse, "\n");
        gf4_poly_t div = gf4_poly_init_zero(10);
        gf4_poly_t rem = gf4_poly_init_zero(10);
        gf4_poly_t m = gf4_poly_init_zero(20);
        gf4_poly_mul(&m, &maybe_inverse, poly);
        gf4_poly_div_rem(&div, &rem, &m, &modulus);
        printf("inv*poly %% modulus: ");
        gf4_poly_pretty_print(&rem, "\n");
        gf4_poly_deinit(&div);
        gf4_poly_deinit(&rem);
        gf4_poly_deinit(&m);
    } else {
        printf("inverse does not exist!\n");
    }

    gf4_poly_deinit(&poly1);
    gf4_poly_deinit(&poly2);
    gf4_poly_deinit(&modulus);
    gf4_poly_deinit(&maybe_inverse);
	return 0;
    */

    encoding_context_t ec;
    decoding_context_t dc;
    // size_t block_size = 2339;
    // size_t block_weight = 37;
    size_t block_size = 7;
    size_t block_weight = 3;
    contexts_init(&ec, &dc, block_size, block_weight);

    gf4_poly_t msg = gf4_poly_init_zero(block_size);
    gf4_poly_t encoded = gf4_poly_init_zero(2*block_size);
    gf4_poly_t syndrome = gf4_poly_init_zero(block_size);

    fprintf(stderr, "sec: ");
    gf4_poly_coeff_print(&ec.second_block_G, block_size, stderr, "\n");
    fprintf(stderr, "h0:  ");
    gf4_poly_coeff_print(&dc.h0, block_size, stderr, "\n");
    fprintf(stderr, "h1:  ");
    gf4_poly_coeff_print(&dc.h1, block_size, stderr, "\n");

    random_weighted_gf4_poly(&msg, block_size, block_weight);
    fprintf(stderr, "msg: ");
    gf4_poly_coeff_print(&msg, block_size, stderr, "\n");
    enc_encode(&encoded, &msg, &ec);
    fprintf(stderr, "enc: ");
    gf4_poly_coeff_print(&encoded, 2*block_size, stderr, "\n");
    dec_calculate_syndrome(&syndrome, &encoded, &dc);
    fprintf(stderr, "syn: ");
    gf4_poly_coeff_print(&syndrome, block_size, stderr, "\n");

    if (gf4_poly_is_zero(&syndrome)) {
        printf("ok!\n");
    } else {
        printf("not ok!\n");
    }

    gf4_poly_deinit(&msg);
    gf4_poly_deinit(&encoded);
    gf4_poly_deinit(&syndrome);
    contexts_deinit(&ec, &dc);
    return 0;
}
