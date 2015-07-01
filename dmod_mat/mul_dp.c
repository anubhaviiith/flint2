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

void _dmod_mat_mul_dp(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B)
{
    #if HAVE_BLAS
    
    slong i, j, z;
    slong m, n, k;

    m = C->nrows;
    n = C->ncols;
    k = A->ncols;

    if (A->ncols != B->nrows)
        return;
   
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            dmod_mat_entry(C, i, j) = _dmod_vec_dot_ld(dmod_mat_entry_ptr(A, i, 0), dmod_mat_entry_ptr(B, 0, j), B->ld,  k, C->mod);
        }
    } 

    #endif
}