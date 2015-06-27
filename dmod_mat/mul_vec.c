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
#include <cblas.h>
#include "dmod_mat.h"

void _dmod_mat_mul_vec(double *res, dmod_mat_t A, const double *x, slong lenx, dmod_t mod)
{
    #if HAVE_BLAS
    slong m, n, i;

    m = A->nrows;
    n = A->ncols;
    
    if (A->ncols != lenx)
        return;
    
    cblas_dgemv(101, 111, m, n, 1, A->rows, n, x, 1, 0, res, 1); 
    
    for (i = 0; i < m; i++)
    {
        res[i] = dmod_reduce(res[i], mod);
    }

    #endif
}
