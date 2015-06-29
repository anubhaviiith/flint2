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

void _dmod_mat_window_init(dmod_mat_t window, const dmod_mat_t A, slong r1, slong c1, slong m, slong n)
{
    #if HAVE_BLAS
    slong i;

    window->rows = A->rows + MATRIX_IDX(A->ld, r1, c1);
    window->nrows = m;
    window->ncols = n;
    window->ld = A->ld;
    window->mod = A->mod;
    #endif

}
