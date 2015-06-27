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

void _dmod_mat_add(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B)
{
    slong i, j;
    
    if ((A->nrows != B->nrows) || (A->ncols != B->ncols))
        return;

    for (i = 0; i < C->nrows; i++)    
    {
        for (j = 0; j < C->ncols; j++)
        {
            C->rows[MATRIX_IDX(C->ncols, i, j)] = 
                dmod_add(A->rows[MATRIX_IDX(A->ncols, i, j)], B->rows[MATRIX_IDX(B->ncols, i, j)], C->mod);
        }
    }
}
