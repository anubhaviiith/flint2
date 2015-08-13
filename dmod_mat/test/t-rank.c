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
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "dmod_vec.h"
#include "nmod_mat.h"
#include "dmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    nmod_mat_t A;
    dmod_mat_t A_d;
    slong i, m, n, d, r, q, w, rand, rank1, rank2, bits;
    FLINT_TEST_INIT(state);
    
    dmod_t mod;

    slong limit_dbl = (1UL << (FLINT_D_BITS/2));
    flint_printf("rank....");
    fflush(stdout);

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
         
        bits = n_randint(state, 27);
        if (bits < 2)
            bits = 2;
        rand = n_randprime(state, bits, 0);
        dmod_init(&mod, rand); 
        

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            nmod_mat_init(A, m, n, rand);
            _dmod_mat_init(A_d, m, n, mod);
            nmod_mat_randrank(A, state, r);
            
            for (q = 0; q < m; q++)
            {
                for (w = 0; w < n; w++)
                {
                    dmod_mat_entry(A_d, q, w) = (double)A->rows[q][w];
                }
            }
            rank1 = nmod_mat_rank(A);
            rank2 = _dmod_mat_rank(A_d);

            if (r != rank2 || rank1 != rank2)
            {
                flint_printf("FAIL:\n");
                gmp_printf("%Mu %Mu\n", rank1, rank2);
                flint_printf("wrong rank!\n");
                abort();
            }
            nmod_mat_clear(A);
            _dmod_mat_clear(A_d);
        }
    }
    
    /* Dense */

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        
        bits = n_randint(state, 27);
        if (bits < 2)
            bits = 2;
        rand = n_randprime(state, bits, 0);
        dmod_init(&mod, rand); 
        
        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            d = n_randint(state, 2*m*n + 1);
            nmod_mat_init(A, m, n, rand);
            _dmod_mat_init(A_d, m, n, mod);
            nmod_mat_randrank(A, state, r);
            
            nmod_mat_randops(A, d, state);
            
            for (q = 0; q < m; q++)
            {
                for (w = 0; w < n; w++)
                {
                    dmod_mat_entry(A_d, q, w) = (double)A->rows[q][w];
                }
            }
            
            rank1 = nmod_mat_rank(A);
            rank2 = _dmod_mat_rank(A_d);
    
            if (rank1 != rank2 || r != rank2)
            {
                gmp_printf("%Mu %Mu %Mu\n", rank1, r, rank2);
                flint_printf("wrong rank!\n");
                abort();
            }
            nmod_mat_clear(A);
            _dmod_mat_clear(A_d);
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
