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


#include <cblas.h>
#include "dmod_vec.h"
#include <stdio.h>
#include "fmpz.h"


double _dmod_vec_dot(const double *vec1, const double *vec2, slong N, dmod_t mod)
{
    double sum = 0;
    slong i, ptr = 0;
    double *a, *b;

    ulong window = n_powmod_precomp(2,(FLINT_D_BITS - (2 * mod.b)), mod.n, mod.ninv);
    
    a = _dmod_vec_init(N);
    b = _dmod_vec_init(N);

    for (i = 0; i < N; i++)
    {
        a[ptr] = dmod_precomp(vec1[i], mod.n, mod.ninv);
        b[ptr] = dmod_precomp(vec2[i], mod.n, mod.ninv);
        ptr++;
        if (i % window == 0)
        {
            sum = dmod_precomp(sum, mod.n, mod.ninv);
            sum += cblas_ddot(ptr ,a, 1, b, 1);
            ptr = 0;
        }
    }

    sum = dmod_precomp(sum, mod.n, mod.ninv);
    sum += cblas_ddot(ptr ,a, 1, b, 1);
    
    _dmod_vec_clear(a);
    _dmod_vec_clear(b);

    return sum;
}
