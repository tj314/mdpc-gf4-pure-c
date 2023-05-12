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

#include "gf4.h"

// helper function
bool gf4_is_in_range(gf4_t a) {
    return a <= GF4_MAX_VALUE;
}

const char * gf4_to_str(gf4_t a) {
    assert(gf4_is_in_range(a));
    switch (a) {
        case 0:
            return "0";
        case 1:
            return "1";
        case 2:
            return "a";
        case 3:
            return "(a+1)";
        default:
            return "";
    }
}


// operations
gf4_t gf4_add(gf4_t a, gf4_t b) {
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
    return a ^ b;
}

gf4_t gf4_mul(gf4_t a, gf4_t b) {
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
    return GF4_MULTIPLICATION[a][b];
}

gf4_t gf4_div(gf4_t a, gf4_t b) {
    assert(gf4_is_in_range(a));
	assert(gf4_is_in_range(b));
	assert(0 != b);
    return GF4_DIVISION[a][b-1];
}