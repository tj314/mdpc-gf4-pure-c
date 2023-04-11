#ifdef RUNTESTS
#include "tests.h"

// TESTS
void test_dec_calculate_syndrome() {
    // setup
    decoding_context_t dc;
    dc.block_size = 3;
    dc.h0 = gf4_poly_init_zero(dc.block_size);
    dc.h1 = gf4_poly_init_zero(dc.block_size);
    gf4_poly_set_coefficient(&dc.h0, 0, 1);
    gf4_poly_set_coefficient(&dc.h1, 0, 2);
    gf4_poly_set_coefficient(&dc.h1, 2, 3);
    gf4_poly_t syndrome = gf4_poly_init_zero(dc.block_size);

    // test 1
    gf4_poly_t vec = gf4_poly_init_zero(2*dc.block_size);
    dec_calculate_syndrome(&syndrome, &vec, &dc);
    assert(0 == syndrome.coefficients[0]);
    assert(0 == syndrome.coefficients[1]);
    assert(0 == syndrome.coefficients[2]);
    assert(0 == vec.coefficients[0]);
    assert(0 == vec.coefficients[1]);
    assert(0 == vec.coefficients[2]);
    assert(0 == vec.coefficients[3]);
    assert(0 == vec.coefficients[4]);
    assert(0 == vec.coefficients[5]);
    assert(1 == dc.h0.coefficients[0]);
    assert(0 == dc.h0.coefficients[1]);
    assert(0 == dc.h0.coefficients[2]);
    assert(2 == dc.h1.coefficients[0]);
    assert(0 == dc.h1.coefficients[1]);
    assert(3 == dc.h1.coefficients[2]);
    assert(3 == dc.block_size);
    assert(0 == dc.h0.degree);
    assert(2 == dc.h1.degree);

    // test 2
    gf4_poly_set_coefficient(&vec, 0, 1);
    gf4_poly_set_coefficient(&vec, 1, 2);
    gf4_poly_set_coefficient(&vec, 4, 1);
    gf4_poly_set_coefficient(&vec, 5, 3);

    dec_calculate_syndrome(&syndrome, &vec, &dc);

    assert(3 == syndrome.coefficients[0]);
    assert(0 == syndrome.coefficients[1]);
    assert(2 == syndrome.coefficients[2]);
    assert(1 == vec.coefficients[0]);
    assert(2 == vec.coefficients[1]);
    assert(0 == vec.coefficients[2]);
    assert(0 == vec.coefficients[3]);
    assert(1 == vec.coefficients[4]);
    assert(3 == vec.coefficients[5]);
    assert(1 == dc.h0.coefficients[0]);
    assert(0 == dc.h0.coefficients[1]);
    assert(0 == dc.h0.coefficients[2]);
    assert(2 == dc.h1.coefficients[0]);
    assert(0 == dc.h1.coefficients[1]);
    assert(3 == dc.h1.coefficients[2]);
    assert(3 == dc.block_size);
    assert(0 == dc.h0.degree);
    assert(2 == dc.h1.degree);

    // cleanup
    gf4_poly_deinit(&vec);
    gf4_poly_deinit(&dc.h0);
    gf4_poly_deinit(&dc.h1);
    gf4_poly_deinit(&syndrome);
}
#endif