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

void _dmod_mat_solve_triu_recursive(dmod_mat_t X, const dmod_mat_t U, const dmod_mat_t B, int unit)
{
    dmod_mat_t UA, UB, UD, XX, XY, BX, BY, tmp;
    slong r, n, m;

    n = U->nrows;
    m = B->ncols;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    _dmod_mat_window_init(UA, U, 0, 0, r, r);
    _dmod_mat_window_init(UB, U, 0, r, r, n - r);
    _dmod_mat_window_init(UD, U, r, r, n - r, n - r);
    _dmod_mat_window_init(BX, B, 0, 0, r, m);
    _dmod_mat_window_init(BY, B, r, 0, n - r, m);
    _dmod_mat_window_init(XX, X, 0, 0, r, m);
    _dmod_mat_window_init(XY, X, r, 0, n - r, m);

    _dmod_mat_solve_triu(XY, UD, BY, unit);
    
    _dmod_mat_init(tmp, r, m, X->mod);
    _dmod_mat_mul(tmp, UB, XY);
    _dmod_mat_sub(XX, BX, tmp);

    _dmod_mat_solve_triu(XX, UA, XX, unit);

    _dmod_mat_window_clear(UA);
    _dmod_mat_window_clear(UB);
    _dmod_mat_window_clear(UD);
    _dmod_mat_window_clear(BX);
    _dmod_mat_window_clear(BY);
    _dmod_mat_window_clear(XX);
    _dmod_mat_window_clear(XY);
}
