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

void _dmod_mat_scalar_mul(dmod_mat_t B, const dmod_mat_t A, double c)
{
    #if HAVE_BLAS
    slong i, j;

    for (i = 0; i < A->nrows; i++)
    {
        for (j = 0; j < A->ncols; j++)
        {
            _dmod_mat_set(B, i, j, dmod_mul(dmod_mat_entry(A, i, j), c, A->mod));
        }
    }
    #endif
}
