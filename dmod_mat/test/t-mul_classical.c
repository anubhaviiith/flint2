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
#include "dmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, i, j, k, rep, rand;
    FLINT_TEST_INIT(state);
    

    flint_printf("mul ....");
    fflush(stdout);

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        dmod_mat_t A, B, C;
        
        m = n_randint(state, 50);
        n = n_randint(state, 50);
        k = n_randint(state, 50);
        rand = n_randint(state, 50);
        while (m==0 || n==0 || k ==0 || rand == 0)
        {
            m = n_randint(state, 50);
            n = n_randint(state, 50);
            k = n_randint(state, 50);
            rand = n_randint(state, 50);

        }
        dmod_t mod;

        dmod_init(&mod, rand); 

        _dmod_mat_init(C, m, n, mod);
        _dmod_mat_init(A, m, k, mod);
        _dmod_mat_init(B, k, n, mod);
        
        _dmod_mat_randtest(A, state);
        _dmod_mat_randtest(B, state);
        _dmod_mat_mul(C, A, B); 
        
        _dmod_mat_clear(A);
        _dmod_mat_clear(B);
        _dmod_mat_clear(C);
        break;
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
