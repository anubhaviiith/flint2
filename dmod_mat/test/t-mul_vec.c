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
    slong m, n, i, j, rep, rand;
    FLINT_TEST_INIT(state);
    

    flint_printf("mul_vec ....");
    fflush(stdout);

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, result;
        dmod_mat_t A_d;
        double *x, *y;
        ulong limit_dbl;

        m = n_randint(state, 100);
        n = n_randint(state, 100);
        
        while (m == 0 || n == 0)
        {
            m = n_randint(state, 100);
            n = n_randint(state, 100);
        }

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
        
        x = _dmod_vec_init(n);
        y = _dmod_vec_init(m);
        
        _dmod_mat_init(A_d, m, n, mod);
       
        nmod_mat_init(A, m, n, modn.n);
        nmod_mat_init(B, n, 1, modn.n);
        nmod_mat_init(result, m, 1, modn.n);
        
        
        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);
        
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                dmod_mat_entry(A_d, i, j) = (double)A->rows[i][j];
            }
        }

        for (i = 0; i < n; i++)
            x[i] = (double)B->rows[i][0];
        
        
        nmod_mat_mul(result, A, B); 
        _dmod_mat_mul_vec(y, A_d, x, n, mod);
        
        for (i = 0; i < m; i++)
        { 
            if (y[i] != (double)result->rows[i][0])
            {
                flint_printf("FAIL\n");
                abort();
            }
            
        }
        nmod_mat_clear(B);
        nmod_mat_clear(result);
        nmod_mat_clear(A);
        
        _dmod_mat_clear(A_d);
        _dmod_vec_clear(x);
        _dmod_vec_clear(y);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
    #endif
}
