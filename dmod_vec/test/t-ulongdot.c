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
#include <gmp.h>

#include "flint.h"
#include "dmod_vec.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

#include "fmpz.h"
#include "nmod_vec.h"

int main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("dot....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        mp_limb_t result1, result2;

        mp_ptr c, d;
        
        fmpz_t m, window, base, limit;

        fmpz_init(m);
        fmpz_init(window);
        fmpz_init(limit);
        fmpz_init(base);
        
        fmpz_set_ui(base, 2);
        
        fmpz_pow_ui(limit, base, FLINT_D_BITS/2);
        ulong limit_ulong = fmpz_get_ui(limit);
        mp_limb_t m_d = n_randint(state, limit_ulong);
        

        dmod_t mod;
        dmod_init(&mod, m_d);
        
        fmpz_pow_ui(window, base, FLINT_D_BITS - 2*mod.b);
        ulong win = fmpz_get_ui(window);
        
        ulong len = 1000;
        
        if (!len)
            continue;

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);

        c = _nmod_vec_init(len);
        d = _nmod_vec_init(len);
        

        for (j = 0; j < len; j++)
        {
            fmpz_randtest_not_zero(a + j, state, FLINT_BITS);
            fmpz_mod_ui(a + j, a + j, mod.n);
            c[j] = fmpz_get_ui(a + j);

        }
        for (j = 0; j < len; j++)
        {
            fmpz_randtest_not_zero(b + j, state, FLINT_BITS);
            fmpz_mod_ui(b + j, b + j, mod.n);
            d[j] = fmpz_get_ui(b + j);
        }
        
        result1 = _dmod_vec_dot(c, d, len, win, mod); 
        
        fmpz_t sum;
        fmpz_init(sum);
        _fmpz_vec_dot(sum, a, b, len);
        fmpz_mod_ui(sum, sum, mod.n);
        result2 = fmpz_get_ui(sum); 
        fmpz_clear(sum);
        
        if(result1 != result2)
        {
            printf("FAIL");
            abort();
        }
        
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    
        fmpz_clear(window);  
        fmpz_clear(m);
        fmpz_clear(limit);
        fmpz_clear(base);

    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
