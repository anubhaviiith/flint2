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
#include "perm.h"

slong _dmod_mat_rref(dmod_mat_t A, slong * pivots_nonpivots, slong * P)
{
    slong i, j, k, n, rank;
    slong * pivots;
    slong * nonpivots;

    dmod_mat_t U, V;

    n = A->ncols;

    rank = _dmod_mat_lu(P, A, 0);

    if (rank == 0)
    {
        for (i = 0; i < n; i++)
            pivots_nonpivots[i] = i;
        return rank;
    }

    /* Clear L */
    for (i = 0; i < A->nrows; i++)
        for (j = 0; j < FLINT_MIN(i, rank); j++)
            dmod_mat_entry(A, i, j) = 0;

    /* We now reorder U to proper upper triangular form U | V
       with U full-rank triangular, set V = U^(-1) V, and then
       put the column back in the original order.

       An improvement for some matrices would be to compress V by
       discarding columns containing nothing but zeros. */

    _dmod_mat_init(U, rank, rank, A->mod);
    _dmod_mat_init(V, rank, n - rank, A->mod);

    pivots = pivots_nonpivots;
    nonpivots = pivots_nonpivots + rank;

    for (i = j = k = 0; i < rank; i++)
    {
        while (dmod_mat_entry(A, i, j) == 0)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }
        pivots[i] = j;
        j++;
    }
    while (k < n - rank)
    {
        nonpivots[k] = j;
        k++;
        j++;
    }

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
            dmod_mat_entry(U, j, i) = dmod_mat_entry(A, j, pivots[i]);
    }

    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
            dmod_mat_entry(V, j, i) = dmod_mat_entry(A, j, nonpivots[i]);
    }

    _dmod_mat_solve_triu(V, U, V, 0);

    /* Clear pivot columns */
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
            dmod_mat_entry(A, j, pivots[i]) = (i == j);
    }

    /* Write back the actual content */
    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
            dmod_mat_entry(A, j, nonpivots[i]) = dmod_mat_entry(V, j, i);
    }

    _dmod_mat_clear(U);
    _dmod_mat_clear(V);

    return rank;
}

slong dmod_mat_rref(dmod_mat_t A)
{
    slong rank, * pivots_nonpivots, * P;
    pivots_nonpivots = flint_malloc(sizeof(slong) * A->ncols);
    P = _perm_init(A->nrows);

    rank = _dmod_mat_rref(A, pivots_nonpivots, P);

    flint_free(pivots_nonpivots);
    _perm_clear(P);

    return rank;
}
