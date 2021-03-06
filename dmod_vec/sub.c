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
#include<stdio.h>

void _dmod_vec_sub(double *result, const double *vec1, const double *vec2, slong len, dmod_t mod)
{
    #if HAVE_BLAS
    slong i;
    
    for (i = 0; i < len; i++)
    {
        result[i] = dmod_sub(vec1[i], vec2[i], mod);
    }
    #endif

}
