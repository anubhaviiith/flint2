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


static __inline__ int _dmod_mat_pivot(dmod_mat_t A, slong * P, slong start_row, slong col)
{
    slong j, t;

    if (dmod_mat_entry(A, start_row, col) != 0)
        return 1;

    for (j = start_row + 1; j < A->nrows; j++)
    {
        if (dmod_mat_entry(A, j, col) != 0)
        {
            _dmod_vec_swap(dmod_mat_entry_ptr(A, j, 0), dmod_mat_entry_ptr(A, start_row, 0), A->ncols);
            t = P[j];
            P[j] = P[start_row];
            P[start_row] = t;

            return -1;
        }
    }
    return 0;
}


slong _dmod_mat_lu_classical(slong * P, dmod_mat_t A, int rank_check)
{

    double d, e;
    
    double *a;

    slong i, m, n, rank, length, row, col;

    m = A->nrows;
    n = A->ncols;

    a = A->rows;/*dmod_mat_entry_ptr(A, 0, 0);
    */
    rank = row = col = 0;

    for (i = 0; i < m; i++)
        P[i] = i;

    while (row < m && col < n)
    {
        if (_dmod_mat_pivot(A, P, row, col) == 0)
        {
            if (rank_check)
                return 0;
            col++;
            continue;
        }

        rank++;

        d = a[(row * A->ld) + col];
        d = dmod_inv(d, A->mod);
        length = n - col - 1;
        
        for (i = row + 1; i < m; i++)
        {
            e = dmod_mul(a[(i * A->ld) + col], d, A->mod);
            if (length != 0)
                _dmod_vec_scalar_addmul_dmod(a + (i * A->ld) + col + 1,
                    a + (row * A->ld) + col + 1, dmod_neg(e, A->mod), length, A->mod);

            a[i*A->ld + col] = 0;
            a[i*A->ld + rank - 1] = e;
        }
        row++;
        col++;
    }
    
    
    return rank;
}
