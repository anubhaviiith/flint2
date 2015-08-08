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
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "dmod_vec.h"

void _dmod_mat_transpose(dmod_mat_t B, const dmod_mat_t A)
{
    double tmp;

    slong i, j;

    if (B->nrows != A->ncols || B->ncols != A->nrows)
    {
        flint_printf("Exception (dmod_mat_transpose). Incompatible dimensions.\n");
        abort();
    }

    if (A == B)
    {
        for (i = 0; i < A->nrows - 1; i++)
        {
            for (j = i + 1; j < A->ncols; j++)
            {
                tmp = dmod_mat_entry(A, i, j);
                dmod_mat_entry(A, i, j) = dmod_mat_entry(A, j, i);
                dmod_mat_entry(A, j, i) = tmp;
            }
        }
    }
    else 
    {
        for (i = 0; i < B->nrows; i++)
            for (j = 0; j < B->ncols; j++)
                dmod_mat_entry(B, i, j) = dmod_mat_entry(A, j, i);
    }
}
