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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"

slong _dmod_mat_nullspace(dmod_mat_t X, const dmod_mat_t A)
{
    slong i, j, k, m, n, rank, nullity;
    slong * p;
    slong * pivots;
    slong * nonpivots;
    dmod_mat_t tmp;

    m = A->nrows;
    n = A->ncols;

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    _dmod_mat_init(tmp, m, n, A->mod);
    _dmod_mat_copy(tmp, A);

    rank = _dmod_mat_rref(tmp);
    nullity = n - rank;

    _dmod_mat_zero(X);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            dmod_mat_entry(X, i, i) = 1;
    }
    else if (nullity)
    {
        pivots = p;            /* length = rank */
        nonpivots = p + rank;  /* length = nullity */

        for (i = j = k = 0; i < rank; i++)
        {
            while (dmod_mat_entry(tmp, i, j) == 0)
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < nullity)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
            {
                double c = dmod_mat_entry(tmp, j, nonpivots[i]);
                dmod_mat_entry(X, pivots[j], i) = dmod_neg(c, A->mod);
            }

            dmod_mat_entry(X, nonpivots[i], i) = 1;
        }
    }

    flint_free(p);
    _dmod_mat_clear(tmp);

    return nullity;
}
