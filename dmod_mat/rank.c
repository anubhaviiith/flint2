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
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include "dmod_mat.h"


slong _dmod_mat_rank(const dmod_mat_t A)
{
    slong m, n, rank;
    slong * perm;
    dmod_mat_t tmp;

    m = A->nrows;
    n = A->ncols;

    if (m == 0 || n == 0)
        return 0;

    _dmod_mat_init(tmp, m, n, A->mod);
    _dmod_mat_copy(tmp, A);

    perm = flint_malloc(sizeof(slong) * m);

    rank = _dmod_mat_lu_classical(perm, tmp, 0);

    flint_free(perm);
    _dmod_mat_clear(tmp);
    return rank;
}
