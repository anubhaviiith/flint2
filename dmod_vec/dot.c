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


double _dmod_vec_dot(const double *vec1, const double *vec2, slong N, ulong window, dmod_t mod)
{
    double sum = 0;
    slong i, j, ptr = 0;
    double *a, *b;
    mp_limb_t res = 0;

    a = _dmod_vec_init(N);
    b = _dmod_vec_init(N);
    printf("%lld\n", window); 
    for (i = 0; i < N; i++)
    {
        a[ptr] = vec1[i];
        b[ptr] = vec2[i];
        ptr++;
        if (i % window == 0)
        {
            sum += cblas_ddot(ptr ,a , 1, b, 1);
            sum = n_mod2_precomp(sum, mod.n, mod.ninv);
            res = n_mod2_precomp(res, mod.n, mod.ninv);
            printf("ptr = %lld i = %lld res = %lld sum = %lf\n", ptr, i, res, sum);
            for (j = 0; j < ptr; ++j)
            {
                printf("%lld ", a[j]);
            }
            printf("\n");
            for (j = 0; j < ptr; ++j)
            {
                printf("%lld ", b[j]);
            }
            printf("\n");


            ptr = 0;
        }
        res += vec1[i]*vec2[i]; 
    }
    
    printf("%lld %lld %lld\n", res, ptr, N);

    sum = n_mod2_precomp(sum, mod.n, mod.ninv);
    sum += cblas_ddot(ptr ,a, 1, b, 1);
    
    _dmod_vec_clear(a);
    _dmod_vec_clear(b);

    return sum;
}
