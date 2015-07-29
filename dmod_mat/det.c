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
#include "perm.h"


double _dmod_mat_det(dmod_mat_t A)
{
    double det;
    slong * P;

    slong m = A->nrows;
    slong rank;
    slong i, j;

    P = flint_malloc(sizeof(slong) * m);
    
    rank = _dmod_mat_lu_classical(P, A, 1);
    
    det = 0;

    if (rank == m)
    {
        det = 1;
        for (i = 0; i < m; i++)
        {
            det = dmod_mul(det, dmod_mat_entry(A, i, i), A->mod);
        }
    }

    if (_perm_parity(P, m) == 1)
        det = dmod_neg(det, A->mod);
    
    flint_free(P);
    return det;
}

mp_limb_t dmod_mat_det(const dmod_mat_t A)
{
    dmod_mat_t tmp;
    double det;
    slong dim = A->nrows;

    if (dim != A->ncols)
    {
        flint_printf("Exception (nmod_mat_det). Non-square matrix.\n");
        abort();
    }

    if (dim == 0) return 1;
    if (dim == 1) return dmod_mat_entry(A, 0, 0);

    _dmod_mat_init(tmp, A->nrows, A->ncols, A->mod);
    _dmod_mat_copy(tmp, A);
    det = _dmod_mat_det(tmp);
    _dmod_mat_clear(tmp);

    return det;
}
