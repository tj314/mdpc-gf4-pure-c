#include <stdio.h>
#include "src/gf4.h"
#include "src/gf4_poly.h"


int main(void) {
	gf4_t a = 2;
	gf4_t b = 1;
	gf4_t c = gf4_add(a, b);
	printf("%s + %s = %s\n", gf4_to_str(a), gf4_to_str(b), gf4_to_str(c));
	return 0;
}
