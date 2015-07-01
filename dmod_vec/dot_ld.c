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

    Copyright (C) 2015 Anubhav Srivastava

******************************************************************************/


#include <cblas.h>
#include "dmod_vec.h"
#include "nmod_vec.h"
#include <stdio.h>
#include "fmpz.h"
#include "ulong_extras.h"
#include <math.h>

double _dmod_vec_dot_ld(const double *vec1, const double *vec2, slong ld, slong len, dmod_t mod)
{
    double res1 = 0.0;
    #if HAVE_BLAS
    slong i;
    double val = 0;
    slong window = (1UL << (FLINT_D_BITS - 2*mod.nbits));
    
    for (i = 0; i < (len - window); i += window)
    {
        val = cblas_ddot(window, vec1 + i, 1, vec2 + i*ld, ld);

        val = dmod_reduce(val, mod);

        res1 = dmod_add(res1, val, mod);
    }
    if (i < len)
    {
        val = cblas_ddot(len - i, vec1 + i, 1, vec2 + i*ld, ld);  
        val = dmod_reduce(val, mod);

        res1 = dmod_add(res1, val, mod);
        
    }
    #endif
    return res1;
}
