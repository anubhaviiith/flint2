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
#include "dmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i, j, k, rep, rand;
    FLINT_TEST_INIT(state);
    

    flint_printf("window init/clear ....");
    fflush(stdout);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        dmod_mat_t A_d, result_d, window;
        ulong limit_dbl;
        mp_limb_t r1, r2, c1, c2, m, n;

        m = n_randint(state, 100);
        n = n_randint(state, 100);
        
        while (m <= 2 || n <= 2)
        {
            m = n_randint(state, 100);
            n = n_randint(state, 100);
        }

        limit_dbl = (1UL << (FLINT_D_BITS - 1));
        rand = n_randint(state, limit_dbl);

        while (rand == 0)
        {
            rand = n_randint(state, limit_dbl);
        }
 
        
        r2 = n_randint(state, m);
        while (r2 <= 1)
            r2 = n_randint(state, m);
        
        r1 = n_randint(state, r2);
        while (r1 == 0)
            r1 = n_randint(state, r2);
        
        c2 = n_randint(state, n);
        while (c2 <= 1)
            c2 = n_randint(state, n);
        
        c1 = n_randint(state, c2);
        while (c1 == 0)
            c1 = n_randint(state, c2);


        dmod_t mod;
        dmod_init(&mod, rand); 
       
        _dmod_mat_init(A_d, m, n, mod);
 
        _dmod_mat_randtest(A_d, state);

        _dmod_mat_window_init(window, A_d, r1, c1, r2, c2);
       
        for (i = 0; i < r2 - r1; i++)
        {
            for (j = 0; j < c2 - c1; j++)
            {
                if (window->entry[i][j] != A_d->entry[i + r1][j + c1])
                {
                    flint_printf("FAIL\n");
                    abort();
                }
            }
        }
        _dmod_mat_window_clear(window);
        _dmod_mat_clear(A_d);
        break;

    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
