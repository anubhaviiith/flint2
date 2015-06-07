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

  Copyright (C) 2014 Abhinav Baid

 ******************************************************************************/

#include "dmod_vec.h"
#include "nmod_vec.h"
#include "fmpz.h"

mp_limb_t _umod_dot_fmpztest(mp_srcptr a, mp_srcptr b, slong len, ulong window, nmod_t mod)
{
    fmpz_t sum, p, q;
    fmpz_init(sum);
    fmpz_init(p);
    fmpz_init(q);
    slong j;

    for (j = 0; j < len; j++)
    {
        fmpz_set_ui(p, a[j]);
        fmpz_set_ui(q, b[j]);
        
        if (j % window == 0)
        {
            fmpz_mod_ui(sum, sum, mod.n);
        }
        fmpz_addmul(sum, p, q);
    }
    mp_limb_t ans = fmpz_get_ui(sum); 
    fmpz_clear(sum);
    fmpz_clear(p);
    fmpz_clear(q);

    return ans;
}
