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
#include "dmod_vec.h"

void _dmod_mat_solve_triu_classical(dmod_mat_t X, const dmod_mat_t U, const dmod_mat_t B, int unit)
{
    slong n, m;
    
    n = U->nrows;
    m = B->ncols;

    if (unit)
    {
        dmod_mat_t temp;
        _dmod_mat_init(temp, B->nrows, B->ncols, B->mod);
        _dmod_mat_copy(temp, B); 
        double alpha = 1;

        cblas_dtrsm(101, 141, 121, 111, 132, n, m, alpha, dmod_mat_entry_ptr(U, 0, 0), U->ld, dmod_mat_entry_ptr(temp, 0, 0), temp->ld);

        _dmod_mat_copy(X, temp);
        _dmod_mat_reduce(X);
        _dmod_mat_clear(temp);
    }
    else
    {
        double *inv, *tmp;

        slong i, j;
        inv = _dmod_vec_init(n);
        for (i = 0; i < n; i++)
            inv[i] = dmod_inv(dmod_mat_entry(U, i, i), U->mod);

        tmp = _dmod_vec_init(n);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
                tmp[j] = dmod_mat_entry(X, j, i);

            for (j = n - 1; j >= 0; j--)
            {
                double s;
                s = _dmod_vec_dot(dmod_mat_entry_ptr(U, j, j + 1), tmp + j + 1, n - j - 1, U->mod);
                s = dmod_sub(dmod_mat_entry(B, j, i), s, U->mod);

                if (!unit)
                    s = dmod_mul(s, inv[j], U->mod);
                tmp[j] = s;
            }

            for (j = 0; j < n; j++)
                dmod_mat_entry(X, j, i) = tmp[j];
        }

        _dmod_vec_clear(tmp);

        _dmod_vec_clear(inv);
    }
}
