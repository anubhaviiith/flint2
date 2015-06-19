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
#include "nmod_vec.h"
#include <stdio.h>
#include "fmpz.h"
#include "ulong_extras.h"
#include <math.h>

double _dmod_vec_dot(const double *vec1, const double *vec2, slong N, dmod_t mod)
{
    #if HAVE_BLAS
    slong i, j;
    double res1 = 0.0, val = 0;
    slong window = pow(2, FLINT_D_BITS - 2*(FLINT_BIT_COUNT(mod.n)));
    
    for (i = 0; i < (N - window); i += window)
    {
        val = cblas_ddot(window, vec1 + i, 1, vec2 + i, 1); 
        val = dmod_reduce(val, mod);

        res1 = dmod_add(res1, val, mod);

    }
    if (i < N)
    {
        val = cblas_ddot(N - i, vec1 + i, 1, vec2 + i, 1);  
        val = dmod_reduce(val, mod);

        res1 = dmod_add(res1, val, mod);
    }
    return res1;
    #endif
}
