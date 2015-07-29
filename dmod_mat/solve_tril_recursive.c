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
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "dmod_vec.h"


void _dmod_mat_solve_tril_recursive(dmod_mat_t X, const dmod_mat_t L, const dmod_mat_t B, int unit)
{
    dmod_mat_t LA, LC, LD, XX, XY, BX, BY, tmp;
    slong r, n, m;

    n = L->nrows;
    m = B->ncols;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    _dmod_mat_window_init(LA, L, 0, 0, r, r);
    _dmod_mat_window_init(LC, L, r, 0, n - r, r);
    _dmod_mat_window_init(LD, L, r, r, n - r, n - r);
    _dmod_mat_window_init(BX, B, 0, 0, r, m);
    _dmod_mat_window_init(BY, B, r, 0, n - r, m);
    _dmod_mat_window_init(XX, X, 0, 0, r, m);
    _dmod_mat_window_init(XY, X, r, 0, n - r, m);

    _dmod_mat_solve_tril(XX, LA, BX, unit);
    _dmod_mat_init(tmp, n - r, m, X->mod);
    _dmod_mat_mul(tmp, LC, XX);
    _dmod_mat_sub(XY, BY, tmp);
    _dmod_mat_clear(tmp);
    _dmod_mat_solve_tril(XY, LD, XY, unit);

    _dmod_mat_window_clear(LA);
    _dmod_mat_window_clear(LC);
    _dmod_mat_window_clear(LD);
    _dmod_mat_window_clear(BX);
    _dmod_mat_window_clear(BY);
    _dmod_mat_window_clear(XX);
    _dmod_mat_window_clear(XY);
}
