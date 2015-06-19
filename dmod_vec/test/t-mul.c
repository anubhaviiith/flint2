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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "dmod_vec.h"
#include "d_vec.h"
#include "ulong_extras.h"

#include "fmpz.h"
#include "fmpz_vec.h"
#include "nmod_vec.h"

int main(void)
{
    #if HAVE_BLAS
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("mulmod....");
    fflush(stdout);

    dmod_t mod;
   
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong m_d, limit_n, a, b;
        double c, d, result1, result2;
        
        limit_n = pow(2, FLINT_D_BITS/2);

        a = n_randint(state, limit_n);
        b = n_randint(state, limit_n);

        m_d = n_randint(state, limit_n);
        
        dmod_init(&mod, m_d);
        
        c = (double)a;
        d = (double)b;

        result1 = n_mulmod_precomp(a, b, mod.n, mod.ninv);
        result2 = dmod_mul(c, d, mod);

        if (result1 != result2)
        {
            printf("FAIL");
            abort();
        }
    }
    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
    #endif
}
