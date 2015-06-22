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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include <float.h>
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include <cblas.h>
#include "dmod_mat.h"

void _dmod_mat_rank_update(dmod_mat_t A, const double *x, const double *y, slong lenx, slong leny, dmod_t mod)
{
    slong m, n;
    m = A->nrows;
    n = A->ncols;
    
    if (A->nrows != lenx || A->ncols != leny)
        return;

    cblas_dger (101, m, n, 1, x, 1, y, 1, A->rows, n); 
}
