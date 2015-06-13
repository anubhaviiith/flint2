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

double _dmod_vec_dot(const double *vec1, const double *vec2, slong N, dmod_t mod)
{
    slong i, j, last_i = 0;
    double res1 = 0.0, val;

    for (i = 0; i < N; i+=mod.window)
    {
        if (i != 0)
        {
            val = cblas_ddot(mod.window, vec1 + i - mod.window, 1, vec2 + i - mod.window, 1); 
            val = dmod_mod_precomp(val, mod);

            res1 += val; 
            res1 = dmod_mod_precomp(res1, mod);
            
            last_i = i;
        }
    }
    
    val = cblas_ddot(N - last_i, vec1 + last_i, 1, vec2 + last_i, 1);  
    val = dmod_mod_precomp(val, mod);
    
    res1 += val;
    res1 = dmod_mod_precomp(res1, mod);
    
    return res1;
}
