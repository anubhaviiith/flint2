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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "dmod_vec.h"

void _dmod_mat_mul_classical(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B)
{
    #if HAVE_BLAS
    
    slong i, j, z;
    slong sA, sB, sC;
    slong m, n, k;

    m = C->nrows;
    n = C->ncols;
    k = A->ncols;

    if (A->ncols != B->nrows)
        return;

    slong limit = (1UL << (FLINT_D_BITS - 2*(A->mod).nbits));

    sA = MATRIX_IDX(A->ld, A->wr1, A->wc1);
    sB = MATRIX_IDX(B->ld, B->wr1, B->wc1);
    sC = MATRIX_IDX(C->ld, C->wr1, C->wc1);

    if (k < limit)
        cblas_dgemm(101, 111, 111, m, n, k, 1.0, A->rows + sA, A->ld, B->rows + sB, B->ld, 0.0, C->rows + sC, C->ld);

    else
    {
        double *temp;
        temp = flint_calloc(k, sizeof(double)); 
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                for (z = 0; z < k; z++)
                {
                    temp[z] = B->rows[sB + MATRIX_IDX(B->ld, z, j)];
                }
                C->rows[sC + MATRIX_IDX(C->ld, i, j)] = _dmod_vec_dot(A->rows + sA + MATRIX_IDX(A->ld, i, 0), temp, k, C->mod);
            }
        }
        flint_free(temp);
    }
    for (i = 0; i < m; i++)    
    {
        for (j = 0; j < n; j++)
            C->rows[sC + MATRIX_IDX(C->ld, i, j)] = dmod_reduce(C->rows[sC + MATRIX_IDX(C->ld, i, j)], C->mod);
    }

    #endif
}
