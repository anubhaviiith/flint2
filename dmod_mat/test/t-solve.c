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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i, j, z;
    FLINT_TEST_INIT(state);
    

    flint_printf("solve....");
    fflush(stdout);

    for (z = 0; z < 100 * flint_test_multiplier(); z++)
    {
        nmod_mat_t A, B, Y;
        dmod_mat_t A_d, B_d, Y_d;
        
        slong row, col;

         
        ulong bits = n_randint(state, 27);
        if (bits < 2)
            bits = 2;
        slong m = n_randprime(state, bits, 0);
        dmod_t mod;
        dmod_init(&mod, m);  

        row = n_randint(state, 100);
        col = n_randint(state, 100);
        
        nmod_mat_init(A, row, row, m);
        nmod_mat_init(B, row, col, m);
        nmod_mat_init(Y, row, col, m);
    
        _dmod_mat_init(A_d, row, row, mod);
        _dmod_mat_init(B_d, row, col, mod);
        _dmod_mat_init(Y_d, row, col, mod);


        nmod_mat_randrank(A, state, row);
        nmod_mat_randtest(B, state);
        
        if (n_randint(state, 2))
            nmod_mat_randops(A, 1+n_randint(state, 1+row*row), state);

        for (i = 0; i < row; i++)
        {
            for (j = 0; j < row; j++)
            {
                _dmod_mat_set(A_d, i, j, (double)A->rows[i][j]);
            }
        }

        for (i = 0; i < row; i++)
        {
            for (j = 0; j < col; j++)
            {
                _dmod_mat_set(B_d, i, j, (double)B->rows[i][j]);
            }
        }

        nmod_mat_solve(Y, A, B);
        
        _dmod_mat_solve(Y_d, A_d, B_d);
        

        for (i = 0; i < row; i++)
        {
            for (j = 0; j < col; j++)
            {
                if (dmod_mat_entry(Y_d, i, j) != (double)Y->rows[i][j])
                {
                    flint_printf("FAIL\n");
                    abort();
                }
            }
        }

         
        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(Y);
        
        _dmod_mat_clear(A_d);
        _dmod_mat_clear(B_d);
        _dmod_mat_clear(Y_d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}