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
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include "dmod_mat.h"


int _dmod_mat_solve(dmod_mat_t X, const dmod_mat_t A, const dmod_mat_t B)
{
    slong i, j, rank, *perm;
    dmod_mat_t LU;
    int result;

    if (A->nrows == 0 || B->ncols == 0)
        return 1;

    _dmod_mat_init(LU, A->nrows, A->ncols, A->mod);
    _dmod_mat_copy(LU, A);

    perm = flint_malloc(sizeof(slong) * A->nrows);
    
    for (i = 0; i < A->nrows; i++)
        perm[i] = i;

    rank = _dmod_mat_lu_classical(perm, LU, 1);

    if (rank == A->nrows)
    {
        dmod_mat_t PB;
        _dmod_mat_window_init(PB, B, 0, 0, B->nrows, B->ncols);

        for (i = 0; i < A->nrows; i++)
        {
            for (j = 0; j < A->ncols; j++)
            {
                dmod_mat_entry(PB, i, j) = dmod_mat_entry(B, perm[i], j);
            }
        }

        _dmod_mat_solve_tril(X, LU, PB, 1);
        _dmod_mat_solve_triu(X, LU, X, 0);

        _dmod_mat_window_clear(PB);
        result = 1;
    }
    else
    {
        result = 0;
    }

    _dmod_mat_clear(LU);
    flint_free(perm);

    return result;
}
