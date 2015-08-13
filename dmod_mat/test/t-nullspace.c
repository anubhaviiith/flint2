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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_vec.h"
#include "dmod_mat.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    slong i, j, z;

    FLINT_TEST_INIT(state);

    flint_printf("nullspace....");
    fflush(stdout);    

    for (z = 0; z < 100 * flint_test_multiplier(); z++)
    {
        dmod_mat_t A, B, ker;
        nmod_mat_t A_n;
        dmod_t mod;
        slong m, n, d, r, nullity, nulrank, rand, bits;

        m = n_randint(state, 30);
        n = n_randint(state, 30);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            d = n_randint(state, 2*m*n + 1);

            bits = n_randint(state, 27);
            if (bits < 2)
                bits = 2;
            rand = n_randprime(state, bits, 0);
            dmod_init(&mod, rand);  

            _dmod_mat_init(A, m, n, mod);
            nmod_mat_init(A_n, m, n, rand);
            _dmod_mat_init(ker, n, n, mod);
            _dmod_mat_init(B, m, n, mod);

            nmod_mat_randrank(A_n, state, r);
            /* Densify */
            if (n_randlimb(state) % 2)
                nmod_mat_randops(A_n, d, state);

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    dmod_mat_entry(A, i, j) = (double)A_n->rows[i][j];
                }
            }
            nullity = _dmod_mat_nullspace(ker, A);
            nulrank = _dmod_mat_rank(ker);

            if (nullity != nulrank)
            {
                flint_printf("FAIL:\n");
                flint_printf("rank(ker) != nullity!\n");
                _dmod_mat_print(A);
                flint_printf("\n");
                abort();
            }

            if (nullity + r != n)
            {
                flint_printf("FAIL:\n");
                flint_printf("nullity + rank != n\n");
                _dmod_mat_print(A);
                flint_printf("\n");
                abort();
            }

            _dmod_mat_mul(B, A, ker);

            if (_dmod_mat_rank(B) != 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("A * ker != 0\n");
                _dmod_mat_print(A);
                flint_printf("\n");
                abort();
            }
            _dmod_mat_clear(A);
            _dmod_mat_clear(ker);
            _dmod_mat_clear(B);
            nmod_mat_clear(A_n);
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
