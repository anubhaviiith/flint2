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

    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/


#include "dmod_mat.h"
#include "dmod_vec.h"

void _dmod_mat_concat_vertical(dmod_mat_t res, const dmod_mat_t mat1, const dmod_mat_t mat2)
{
    slong i;
    slong r1 = mat1->nrows;
    slong c1 = mat1->ncols;
    slong r2 = mat2->nrows;

    for (i = 0; i < r1; i++)
        _dmod_vec_copy(dmod_mat_entry_ptr(mat1, i, 0), dmod_mat_entry_ptr(res, i, 0), c1);

    for (i = 0; i < r2; i++)
        _dmod_vec_copy(dmod_mat_entry_ptr(mat2, i, 0), dmod_mat_entry_ptr(res, (i + r1), 0), c1);

}
