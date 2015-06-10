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

mp_limb_t _umod_vec_dot(mp_srcptr vec1, mp_srcptr vec2, slong N, ulong window, umod_t mod)
{
    slong i, j, ptr = 0;
    mp_limb_t res = 0;

    for (i = 0; i < N; i++)
    {
        if (i % window == 0)
        {
            res = dmod_mod_precomp(res, mod);
        }
        res += vec1[i]*vec2[i]; 
    }
    res = dmod_mod_precomp(res, mod);
    return res;
}
