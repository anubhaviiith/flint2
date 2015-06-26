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

    Copyright (C) 2008, 2009 William Hart.
    Copyright (C) 2008, Richard Howell-Peak
    Copyright (C) 2008, Martin Albrecht
    Copyright (C) 2010, Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"

void _dmod_mat_window_init(dmod_mat_t window, const dmod_mat_t mat, slong r1, slong c1, slong r2, slong c2)
{
    slong i, j;
    slong m = r2 -r1, n = c2 -c1; 
    
    
    window->rows = flint_calloc((m * n), sizeof(double));
    window->entry = flint_malloc(m * sizeof(double *));
    
    for (i = 0; i < m; i++)
    {
        window->entry[i] = window->rows + i * n;
    }
    
    for (i = r1; i < r2; i++)
    {
        for (j = c1; j < c2; j++)
        {
            window->entry[i - r1][j - c1] = mat->rows[MATRIX_IDX(mat->ncols, i, j)]; 
        }
    }

    window->nrows = m;
    window->ncols = n;
    window->mod = mat->mod;
}
