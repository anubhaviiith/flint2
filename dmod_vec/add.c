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
#include <stdio.h>
#include <immintrin.h>

void _dmod_vec_add(double *result, const double *vec1, const double *vec2, slong len, dmod_t mod)
{
    #if HAVE_BLAS
    slong i;
    
    for (i = 0; i < len; i++)
    {
        result[i] = dmod_add(vec1[i], vec2[i], mod);
    }
    /*
    for (i = 0; i < (len - 4); i += 4)
    {
        __m256d kA2   = _mm256_load_pd( &vec1[i] );
        __m256d kB2   = _mm256_load_pd( &vec2[i] );

        __m256d kRes = _mm256_add_pd( kA2, kB2 );
        _mm256_stream_pd( &result[i], kRes );
    }
    if (i < len && (i - 3) >= 0)
    {
        result[i - 3] = vec1[i - 3] + vec2[i - 3];
        result[i - 2] = vec1[i - 2] + vec2[i - 2];
        result[i - 1] = vec1[i - 1] + vec2[i - 1];
        result[i] = vec1[i] + vec2[i];
    }
    _dmod_vec_reduce(result, result, len, mod);
    */
#endif
}
