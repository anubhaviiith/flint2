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

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong limit_ulong, m_d;
        fmpz_t a, b, product, base, limit;
        
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(product);
        fmpz_init(base);
        fmpz_init(limit);

        fmpz_randtest_unsigned(a, state, FLINT_BITS);
        fmpz_randtest_unsigned(b, state, FLINT_BITS);
       
        fmpz_set_ui(base, 2);        
        fmpz_pow_ui(limit, base, FLINT_D_BITS/2);
        limit_ulong = fmpz_get_ui(limit);
        m_d = n_randint(state, limit_ulong);
        
        dmod_t mod;
        dmod_init(&mod, m_d);
        
        fmpz_mod_ui(a, a, mod.n);
        fmpz_mod_ui(b, b, mod.n);
    
        fmpz_mul(product, a, b);
        fmpz_mod_ui(product, product, mod.n);

        double result1 = fmpz_get_d(product); 

        double c = fmpz_get_d(a);
        double d = fmpz_get_d(b);

        double result2 = dmod_mulmod_precomp(c, d, mod);
        
        if(result1 != result2)
        {
            printf("FAIL");
            abort();
        }
        
        fmpz_clear(product);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(base);
        fmpz_clear(limit);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
