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
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include "dmod_mat.h"


int _dmod_mat_inv(dmod_mat_t B, const dmod_mat_t A)
{
    dmod_mat_t I;
    slong i, dim;
    int result;

    dim = A->nrows;

    switch (dim)
    {
        case 0:
            result = 1;
            break;

        case 1:
            if (dmod_mat_entry(A, 0, 0) == 0)
            {
                result = 0;
            }
            else
            {
                dmod_mat_entry(B, 0, 0) = dmod_inv(dmod_mat_entry(A, 0, 0), B->mod);
                result = 1;
            }
            break;

        default:
            _dmod_mat_init(I, dim, dim, B->mod);
            for (i = 0; i < dim; i++)
                dmod_mat_entry(I, i, i) = 1;
            result = _dmod_mat_solve(B, A, I);
            _dmod_mat_clear(I);
    }

    return result;
}
