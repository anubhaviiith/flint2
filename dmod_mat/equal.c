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

int _dmod_mat_equal(const dmod_mat_t mat1, const dmod_mat_t mat2)
{
    slong i, j;

    if (mat1->nrows != mat2->nrows || mat1->ncols != mat2->ncols)
        return 0;

    if (mat1->nrows == 0 || mat1->ncols == 0)
        return 1;

    for (i = 0; i < mat1->nrows; i++)
    {
        for (j = 0; j < mat1->ncols; j++)
        {
            if (dmod_mat_entry(mat1, i, j) != dmod_mat_entry(mat2, i, j))
                return 0;
        } 
    }

    return 1;
}
