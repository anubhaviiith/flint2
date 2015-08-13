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
#include "dmod_vec.h"
#include "dmod_mat.h"
#include "ulong_extras.h"
#include "perm.h"

int check_rref_form(slong * perm, dmod_mat_t A, slong rank)
{
    slong i, j, k, prev_pivot;

    /* bottom should be zero */
    for (i = rank; i < A->nrows; i++)
        for (j = 0; j < A->ncols; j++)
            if (dmod_mat_entry(A, i, j) != 0)
                return 0;

    prev_pivot = -1;

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < A->ncols; j++)
        {
            if (dmod_mat_entry(A, i, j) != 0)
            {
                /* pivot should have a higher column index than previous */
                if (j <= prev_pivot)
                    return 0;

                /* column should be 0 ... 0 1 0 ... 0 */
                for (k = 0; k < rank; k++)
                    if (dmod_mat_entry(A, k, j) != (i == k))
                        return 0;

                prev_pivot = j;
                break;
            }
        }
    }

    return 1;
}

int
main(void)
{
    slong i;

    FLINT_TEST_INIT(state);

    flint_printf("rref....");
    fflush(stdout);    

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        dmod_mat_t A, B, C, D;
        dmod_t mod;
        slong j, k, m, n, rank1, rank2, bits, rand;
        slong *perm;
        int equal;
        mp_limb_t c;
 
        bits = n_randint(state, 27);
        if (bits < 2)
            bits = 2;
        rand = n_randprime(state, bits, 0);
        dmod_init(&mod, rand);  

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        perm = _perm_init(2*m);

        _dmod_mat_init(A, m, n, mod);
        _dmod_mat_init(D, 2*m, n, mod);

        _dmod_mat_randtest(A, state);
        _dmod_mat_init(B, m, n, mod);
        _dmod_mat_copy(B, A);
        _dmod_mat_init(C, m, n, mod);
        _dmod_mat_copy(C, A);

        rank1 = dmod_mat_rref(B);

        if (!check_rref_form(perm, B, rank1))
        {
            flint_printf("FAIL (malformed rref)\n");
            _dmod_mat_print(A); flint_printf("\n\n");
            _dmod_mat_print(B); flint_printf("\n\n");
            abort();
        }

        /* Concatenate the original matrix with the rref, scramble the rows,
           and check that the rref is the same */
        _perm_randtest(perm, 2 * m, state);

        for (j = 0; j < m; j++)
        {
            do { c = n_randint(state, mod.n); } while (c == 0);
            for (k = 0; k < n; k++)
                dmod_mat_entry(D, perm[j], k) =
                    dmod_mul(dmod_mat_entry(A, j, k), c, A->mod);
        }

        for (j = 0; j < m; j++)
        {
            do { c = n_randint(state, mod.n); } while (c == 0);
            for (k = 0; k < n; k++)
                dmod_mat_entry(D, perm[m + j], k) =
                    dmod_mul(dmod_mat_entry(B, j, k), c, A->mod);
        }

        rank2 = dmod_mat_rref(D);
        equal = (rank1 == rank2);

        if (equal)
        {
            for (j = 0; j < rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal && (dmod_mat_entry(B, j, k) ==
                                        dmod_mat_entry(D, j, k));
            for (j = rank2; j < 2 * rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal && (dmod_mat_entry(D, j, k) == 0);
        }

        if (!equal)
        {
            flint_printf("FAIL (rank1 = %wd, rank2 = %wd)!\n", rank1, rank2);
            _dmod_mat_print(A); flint_printf("\n\n");
            _dmod_mat_print(B); flint_printf("\n\n");
            _dmod_mat_print(D); flint_printf("\n\n");
            abort();
        }

        _perm_clear(perm);
        _dmod_mat_clear(A);
        _dmod_mat_clear(B);
        _dmod_mat_clear(C);
        _dmod_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
