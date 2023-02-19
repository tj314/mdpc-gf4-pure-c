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
    gf4_poly_t maybe_inverse;

    gf4_poly_set_coefficient(&poly1, 2, 1);
    gf4_poly_set_coefficient(&poly1, 1, 1);
    gf4_poly_set_coefficient(&poly1, 0, 1);

    gf4_poly_set_coefficient(&poly2, 4, 2);
    gf4_poly_set_coefficient(&poly2, 1, 2);

    gf4_poly_set_coefficient(&modulus, 0, 1);
    gf4_poly_set_coefficient(&modulus, 8, 1);

    bool inverse_exists = gf4_poly_invert_slow_byref(&maybe_inverse, &poly1, &modulus);

    gf4_poly_pretty_print(&poly1, "\n");
    gf4_poly_pretty_print(&poly2, "\n");

    if (inverse_exists) {
        printf("inverse: ");
        gf4_poly_pretty_print(&maybe_inverse, "\n");
    } else {
        printf("inverse does not exist!\n");
    }

    gf4_poly_deinit(&poly1);
    gf4_poly_deinit(&poly2);
    gf4_poly_deinit(&modulus);
    gf4_poly_deinit(&maybe_inverse);
	return 0;
}
