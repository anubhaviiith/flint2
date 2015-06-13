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

double _dmod_vec_dot(const double *vec1, const double *vec2, slong N, ulong window, dmod_t mod)
{
    slong i, j, last_i = 0;
    double res1 = 0.0;
    double *v1;
    double *v2;
    v1 = _dmod_vec_init(window);
    v2 = _dmod_vec_init(window);

    for (i = 0; i < N; i+=window)
    {
        if (i != 0 && i % window == 0)
        {
            for (j = i - window; j < i; j++)
            {
                v1[j % window] = vec1[j];
                v2[j % window] = vec2[j];
            } 
            res1 += cblas_ddot(window, v1, 1, v2, 1);  
            last_i = i;

            res1 = dmod_mod_precomp(res1, mod);
        }
    }
    for (j = last_i; j < N; j++)
    {
        v1[j - last_i] = vec1[j];
        v2[j - last_i] = vec2[j];
    } 
    res1 += cblas_ddot(N - last_i, v1, 1, v2, 1);  
    res1 = dmod_mod_precomp(res1, mod);
    return res1;
}
