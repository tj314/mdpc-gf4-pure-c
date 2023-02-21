#include <stdio.h>
#include "src/gf4.h"
#include "src/gf4_poly.h"


int main(void) {
    /*
    PolynomialGF2N<GF4> p1; // invertible
    p1.set_coefficient(0, 1);
    p1.set_coefficient(1, 1);
    p1.set_coefficient(2, 1);

    PolynomialGF2N<GF4> p2; // not invertible
    p2.set_coefficient(1, 2);
    p2.set_coefficient(4, 2);

    PolynomialGF2N<GF4> modulus;
    modulus.set_coefficient(0, 1);
    modulus.set_coefficient(8, 1);
    */



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
}
