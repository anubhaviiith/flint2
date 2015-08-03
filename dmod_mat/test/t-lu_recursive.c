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
    slong i, j;

    FLINT_TEST_INIT(state);
    

    flint_printf("lu_recursive....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, LU;
        dmod_mat_t A_d, LU_d;
        slong rand;
        slong m, n, r, d, rank_dmod, rank_nmod, q, w, limit_dbl;
        slong * P;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        m = 5;
        n = 5;
        
        limit_dbl = (1UL << (FLINT_D_BITS/2));
        rand = n_randint(state, limit_dbl);
        rand = 7;
        dmod_t mod;
        dmod_init(&mod, rand); 
        
        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {

            nmod_mat_init(A, m, n, rand);
            nmod_mat_randrank(A, state, r);
         
            nmod_mat_init_set(LU, A);
            
            _dmod_mat_init(A_d, m, n, mod);
            _dmod_mat_init(LU_d, m, n, mod);

            for (q = 0; q < m; q++)
            {
                for (w = 0; w < n; w++)
                {
                    dmod_mat_entry(A_d, q, w) = (double)A->rows[q][w];
                }
            }
            _dmod_mat_copy(LU_d, A_d);

            
            P = flint_malloc(sizeof(slong) * m);
            rank_dmod = _dmod_mat_lu_recursive(P, LU_d, 0); 
            flint_free(P);
    
            P = flint_malloc(sizeof(slong) * m);
            rank_nmod = nmod_mat_lu_recursive(P, LU, 0);
            flint_free(P);
            
            _dmod_mat_print(LU_d);
            nmod_mat_print_pretty(LU);

            printf("#######################\n");
            if (rank_dmod != rank_nmod || r != rank_dmod)
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                flint_printf("A:");
                _dmod_mat_print(A_d);
                flint_printf("LU:");
                _dmod_mat_print(LU_d);
                abort();
            }
            
            for (q = 0; q < m; q++)
            {
                for (w = 0; w < n; w++)
                {
                    if (dmod_mat_entry(LU_d, q, w) != (double)LU->rows[q][w])
                    {
                        flint_printf("FAIL\n");
                        abort();
                    }
                }
            }

            nmod_mat_clear(A);
            nmod_mat_clear(LU);
            
            _dmod_mat_clear(A_d);
            _dmod_mat_clear(LU_d);
            

        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
