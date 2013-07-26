/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

#include <iostream>
#include <sstream>
#include <string>

#include "fmpz_poly_qxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_manipulation()
{
    fmpz_poly_qxx f;
    f.numref() = "3  0 1 1";
    f.denref() = "3  0 -1 1";
    tassert(!f.is_canonical());
    f.canonicalise();
    tassert(f.numref().to_string() == "2  1 1");
    tassert(f.denref().to_string() == "2  -1 1");
}

void
test_assignment_conversion()
{
    fmpz_poly_qxx f, g;
    f = 1;
    tassert(f.is_one());
    g = 0;
    f = g;
    tassert(f.is_zero());
    tassert(f.to_string() == "0");

    f = "4  1 0 0 1";
    tassert(f.numref().to_string() == "4  1 0 0 1");
    tassert(f.denref().to_string() == "1  1");

    f.denref() = "2  -1 1";
    tassert(f.to_string() == "4  1 0 0 1/2  -1 1");
    g = "4  1 0 0 1/2  -1 1";
    tassert(f == g);

    tassert(f.pretty("x") == "(x^3+1)/(x-1)");
}

void
test_arithmetic()
{
    fmpz_poly_qxx f, g;
    g = "4  1 0 0 1/2  -1 1";
    f = "1  1";

    tassert((f + g).to_string() == "4  0 1 0 1/2  -1 1");
    tassert((g - f).to_string() == "4  2 -1 0 1/2  -1 1");
    tassert(g - f == g + (-f));
    tassert(inv(g).to_string() == "2  -1 1/4  1 0 0 1");

    tassert(2 * g == g * 2);
    f = 2*g;
    tassert(f.numref() == 2*g.numref() && f.denref() == g.denref());
    f /= 2;
    tassert(f == g);

#if 0
    tassert((f * fmpzxx(2)) / fmpzxx(3) == f * fmpqxx(2, 3u));
    tassert((f * fmpzxx(3)) / fmpzxx(2) == f / fmpqxx(2, 3u));
    tassert(f * fmpzxx(2) == f * 2 && f / fmpzxx(2) == f / 2);
#endif

    f = "1  1";
    tassert((f*g).evaluate().numref() == f.numref()*g.numref());
    tassert((f*g).evaluate().denref() == f.denref()*g.denref());
    tassert((f/g).evaluate().numref() == f.numref()*g.denref());
    tassert((f/g).evaluate().denref() == f.denref()*g.numref());
}

// Won't compile if the expression is not done using addmul
template<class T>
bool is_ternary(const T&)
{
    return T::ev_traits_t::temp_rule_t::TERNARY_OP_MARKER + 1;
}

// test stuff which we should get automatically - addmul, references etc
void
test_extras()
{
    // TODO
}

void
test_functions()
{
    fmpz_poly_qxx f;
    tassert(f.is_zero() && !f.is_one());
    f.numref() = 1;
    tassert(f.is_one());

    f = "4  1 0 0 1/2  -1 1";
    tassert(pow(f, 4u) == f*f*f*f);

    tassert(derivative(f).to_string() == "4  -1 0 -3 2/3  1 -2 1");

    // test static methods
    frandxx rand;
    tassert(fmpz_poly_qxx::randtest(rand, 10, 8, 10, 8).numref().degree() < 10);
    tassert(fmpz_poly_qxx::randtest(rand, 10, 8, 10, 8).denref().degree() < 10);
    tassert(flog(height(fmpz_poly_qxx::randtest(
                        rand, 10, 8, 10, 8).numref()), 2) < 8);
    tassert(flog(height(fmpz_poly_qxx::randtest(
                        rand, 10, 8, 10, 8).denref()), 2) < 8);
    tassert(!fmpz_poly_qxx::randtest_not_zero(rand, 10, 8, 10, 8).is_zero());
}

int
main()
{
    std::cout << "fmpz_polyxx....";

    test_manipulation();
    test_assignment_conversion();
    test_arithmetic();
    test_functions();
    test_extras();

    std::cout << "PASS" << std::endl;
    return 0;
}

