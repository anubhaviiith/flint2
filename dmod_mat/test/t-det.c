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

    Copyright (C) 2015 Anubhav Srivastava

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
    #if HAVE_BLAS
    
    slong m, n, i, j, k, rep, rand;
    FLINT_TEST_INIT(state);
    

    flint_printf("det ....");
    fflush(stdout);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A;
        dmod_mat_t A_d;

        slong result1;
        double result2;
        ulong limit_dbl;

        m = n_randint(state, 50);

        limit_dbl = (1UL << (FLINT_D_BITS/2));
        rand = n_randint(state, limit_dbl);
        
        while (rand == 0)
        {
            rand = n_randint(state, limit_dbl);
        }
        dmod_t mod;
        nmod_t modn;
        dmod_init(&mod, rand); 
        nmod_init(&modn, rand); 
        
        _dmod_mat_init(A_d, m, m, mod);
       
        nmod_mat_init(A, m, m, modn.n);
       
        nmod_mat_randtest(A, state);
        
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                _dmod_mat_set(A_d, i, j, (double)A->rows[i][j]);
            }
        }
 
        result1 = nmod_mat_det(A); 
        result2 = dmod_mat_det(A_d);
         

        if((double)result1 != result2)
        {
            printf("FAIL\n");
            abort();
        }

        nmod_mat_clear(A);

        _dmod_mat_clear(A_d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
    #endif
}
