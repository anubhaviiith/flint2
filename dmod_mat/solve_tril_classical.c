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

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "dmod_vec.h"

void _dmod_mat_solve_tril_classical(dmod_mat_t X, const dmod_mat_t L, const dmod_mat_t B, int unit)
{
    slong i, j, n, m;
    
    double *inv, *tmp;

    n = L->nrows;
    m = B->ncols;

    if (!unit)
    {
        inv = _dmod_vec_init(n);
        for (i = 0; i < n; i++)
            inv[i] = dmod_inv(dmod_mat_entry(L, i, i), L->mod);
    }
    else
        inv = NULL;

    tmp = _dmod_vec_init(n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = dmod_mat_entry(X, j, i);

        for (j = 0; j < n; j++)
        {
            double s;
            s = _dmod_vec_dot(dmod_mat_entry_ptr(L, j, 0), tmp, j, L->mod);
            s = dmod_sub(dmod_mat_entry(B, j, i), s, L->mod);

            if (!unit)
                s = dmod_mul(s, inv[j], L->mod);
            tmp[j] = s;
        }

        for (j = 0; j < n; j++)
            dmod_mat_entry(X, j, i) = tmp[j];
    }

    _dmod_vec_clear(tmp);
    
    if (!unit)
        _dmod_vec_clear(inv);
}
