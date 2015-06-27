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

    Copyright (C) 2015 Anubhav Srivastava

******************************************************************************/



#include <gmp.h>
#include <stdlib.h>
#include <float.h>
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include "dmod_mat.h"

void _dmod_mat_mul(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B)
{
    #if HAVE_BLAS
    slong m, n, k, i, j;
    m = A->nrows;
    n = B->ncols;
    k = A->ncols;

    if (A->ncols != B->nrows)
        return;

    cblas_dgemm(101, 111, 111, m, n, k, 1.0, A->rows, k, B->rows, n, 0.0, C->rows, n);
    for (i = 0; i < C->nrows; i++)    
    {
        for (j = 0; j < C->ncols; j++)
            C->rows[MATRIX_IDX(C->ncols, i, j)] = dmod_reduce(C->rows[MATRIX_IDX(C->ncols, i, j)], A->mod);
    }
    #endif
    
}
