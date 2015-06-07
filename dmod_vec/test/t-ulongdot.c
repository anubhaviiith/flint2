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

    flint_printf("dot....");
    fflush(stdout);

    nmod_t mod1;
    nmod_init(&mod1, DBL_MAX);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mp_ptr a, b;
        mp_limb_t m, result1, result2; 
        
        ulong len = n_randint(state, 1000);
        m = n_powmod(2, FLINT_D_BITS/2, mod1.n);
        

        nmod_t mod;
        nmod_init(&mod, m);
        /*gmp_printf("n = %Mu %Mu\n", m, mod.b);
        */
        mp_limb_t window = n_powmod(2, FLINT_D_BITS - 2*mod.b, mod1.n);

            
        if (!len)
            continue;

        a = _nmod_vec_init(len);
        b = _nmod_vec_init(len);

        _umod_vec_randtest(a, state, len, mod);
        _umod_vec_randtest(b, state, len, mod);
 
        /*for( j = 0; j < len; ++j)
        {
            gmp_printf("a = %Mu b = %Mu\n", a[j], b[j]);
        }*/

        /*Code to be tested */
        result1 = _umod_vec_dot(a, b, len, window, mod); /* returns double */

        /* fmpz test */

        result2 = _umod_dot_fmpztest(a, b, len, window, mod);
        gmp_printf("%Mu %Mu\n", result1, result2);
        if(result1 != result2)
        {
            printf("FAIL");
            abort();
        }
        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
