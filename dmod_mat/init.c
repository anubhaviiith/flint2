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

#include "dmod_vec.h"
#include "dmod_mat.h"

void _dmod_mat_init(dmod_mat_t A, slong m, slong n, dmod_t mod)
{
    #if HAVE_BLAS
    slong i;

    A->rows = flint_calloc((m*n), sizeof(double)); 
    A->nrows = m;
    A->ncols = n;
    A->parent = A->rows;
    A->wc1 = 0;
    A->wr1 = 0;
    A->iswin = 0;
    A->ld = n;
    _dmod_mat_set_mod(A, mod.n);
    #endif
}
