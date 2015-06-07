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
#include "fmpz.h"

double _dmod_dot_fmpztest(const double *a, const double *b, slong len, dmod_t mod)
{
    fmpz_t sum, p, q;
    fmpz_init(sum);
    fmpz_init(p);
    fmpz_init(q);
    slong j;

    for (j = 0; j < len; j++)
    {
        fmpz_set_d(p, a[j]);
        fmpz_set_d(q, b[j]);

        if ( j != 0 && j % (n_powmod_precomp(2,(FLINT_D_BITS - (2 * mod.b)), mod.n, mod.ninv)) == 0)
        {
            fmpz_mod_ui(sum, sum, mod.n);
        }
        fmpz_addmul(sum, p, q);
    }
    double ans = fmpz_get_d(sum);
    fmpz_clear(sum);
    fmpz_clear(p);
    fmpz_clear(q);

    return ans;
}
