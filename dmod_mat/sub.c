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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "dmod_vec.h"

void _dmod_mat_sub(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B)
{
    #if HAVE_BLAS
    
    slong i, j;
    
    if (A->iswin == 1 || B->iswin == 1 || C->iswin == 1)
    {
        slong m = C->nrows;
        slong n = C->ncols;

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                C->parent[C->ld * (i + C->wr1) + j + C->wc1] =  
                    dmod_sub(A->parent[A->ld * (i + A->wr1) + j + A->wc1], B->parent[B->ld * (i + B->wr1) +  j + B->wc1], C->mod);
            }
        }
    }
    else
    {
        for (i = 0; i < C->nrows; i++)    
        {
            for (j = 0; j < C->ncols; j++)
            {
                C->rows[MATRIX_IDX(C->ncols, i, j)] = 
                    dmod_sub(A->rows[MATRIX_IDX(A->ncols, i, j)], B->rows[MATRIX_IDX(B->ncols, i, j)], C->mod);
            }
        }

    }
    #endif
}
