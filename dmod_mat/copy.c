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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"

void _dmod_mat_copy(dmod_mat_t B, const dmod_mat_t A)
{
    #if HAVE_BLAS
    slong i;

    if (A->nrows != B->nrows || A->ncols != B->ncols)
        return;

    for (i = 0; i < A->nrows; i++)
    {
        _dmod_vec_copy(dmod_mat_entry_ptr(A, i, 0), dmod_mat_entry_ptr(B, i, 0), A->ncols);
    }
    #endif

}
