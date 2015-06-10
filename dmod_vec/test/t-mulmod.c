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
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("mulmod....");
    fflush(stdout);

    nmod_t mod1;
    nmod_init(&mod1, DBL_MAX);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m, result1, result2; 
        
        mp_limb_t a = n_randint(state, DBL_MAX);
        mp_limb_t b = n_randint(state, DBL_MAX);
        
        m = n_randint(state, n_powmod(2, FLINT_D_BITS/2 - 1, mod1.n));
        
        umod_t mod;
        umod_init(&mod, m);
        

        double c = n_mod2_precomp(a, mod.n, mod.ninv);
        double d = n_mod2_precomp(b, mod.n, mod.ninv);

        
        fmpz_t product, p, q;

        fmpz_init(product);
        fmpz_init(p);
        fmpz_init(q);

        fmpz_set_d(p, c);
        fmpz_set_d(q, d);

        fmpz_mul(product, p, q);

        fmpz_mod_ui(product, product, mod.n);

        result1 = fmpz_get_d(product);  
        
        result2 = dmod_mulmod_precomp(c, d, mod);
        
        /*
        flint_printf("%lld %lld\n", result1, result2);
        */

        if(result1 != result2)
        {
            printf("FAIL");
            abort();
        }
        
        fmpz_clear(product);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
