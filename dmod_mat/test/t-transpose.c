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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, mod, mod2, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("transpose....");
    fflush(stdout);

    /* Rectangular transpose, same modulus */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        dmod_mat_t A, B, C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        ulong limit_dbl = (1UL << (FLINT_D_BITS/2));
        slong rand = n_randint(state, limit_dbl);
        
        dmod_t mod;
        dmod_init(&mod, rand);

        _dmod_mat_init(A, m, n, mod);
        _dmod_mat_init(B, n, m, mod);
        _dmod_mat_init(C, m, n, mod);

        _dmod_mat_randtest(A, state);
        _dmod_mat_randtest(B, state);

        _dmod_mat_transpose(B, A);
        _dmod_mat_transpose(C, B);
        if (!_dmod_mat_equal(C, A))
        {
            flint_printf("FAIL: C != A\n");
            abort();
        }

        _dmod_mat_clear(A);
        _dmod_mat_clear(B);
        _dmod_mat_clear(C);
    }

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        dmod_mat_t A, AT, B, BT, AT2;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        
        ulong limit_dbl = (1UL << (FLINT_D_BITS/2));
        
        slong rand = n_randint(state, limit_dbl);
        
        dmod_t mod;
        dmod_init(&mod, rand);

        slong rand2 = n_randint(state, limit_dbl);
        
        dmod_t mod2;
        dmod_init(&mod2, rand);


        _dmod_mat_init(A, m, n, mod);
        _dmod_mat_init(AT, n, m, mod);
        _dmod_mat_init(B, m, n, mod2);
        _dmod_mat_init(BT, n, m, mod2);
        _dmod_mat_init(AT2, n, m, mod2);

        _dmod_mat_randtest(A, state);
        _dmod_mat_copy(B, A);

        _dmod_mat_transpose(AT, A);
        _dmod_mat_transpose(BT, B);

        _dmod_mat_copy(AT2, AT);

        if (!_dmod_mat_equal(BT, AT2))
        {
            flint_printf("FAIL: AT != BT\n");
            abort();
        }

        _dmod_mat_clear(A);
        _dmod_mat_clear(AT);
        _dmod_mat_clear(AT2);
        _dmod_mat_clear(B);
        _dmod_mat_clear(BT);
    }
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        dmod_mat_t A, B;

        m = n_randint(state, 20);
 
        ulong limit_dbl = (1UL << (FLINT_D_BITS/2));
        slong rand = n_randint(state, limit_dbl);
        
        dmod_t mod;
        dmod_init(&mod, rand);


        _dmod_mat_init(A, m, m, mod);
        _dmod_mat_init(B, m, m, mod);

        _dmod_mat_randtest(A, state);
        _dmod_mat_copy(B, A);

        _dmod_mat_transpose(B, B);
        _dmod_mat_transpose(B, B);

        if (!_dmod_mat_equal(B, A))
        {
            flint_printf("FAIL: B != A\n");
            abort();
        }

        _dmod_mat_clear(A);
        _dmod_mat_clear(B);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
