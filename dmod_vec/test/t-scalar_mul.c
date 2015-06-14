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

    flint_printf("scalar_mul....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz *a, *ans, *result, *result1;
        double *c;

        mp_limb_t limit_ulong, m_d, len, limit_low_ulong;
        fmpz_t m, base, limit, sum, sum1, x, limit_x, mod_fmpz, limit_low;
        
        double alpha;

        fmpz_init(m);
        fmpz_init(x);
        fmpz_init(limit);
        fmpz_init(limit_low);
        fmpz_init(limit_x);
        fmpz_init(base);
        fmpz_init(sum);
        fmpz_init(sum1);
        fmpz_init(mod_fmpz);
        
        fmpz_set_ui(base, 2);        
        fmpz_pow_ui(limit, base, FLINT_D_BITS/2);
        fmpz_pow_ui(limit_low, base, FLINT_D_BITS/4);
        
        limit_ulong = fmpz_get_ui(limit);
        limit_low_ulong = fmpz_get_ui(limit_low);

        if (n_randint(state, 2))
            m_d = n_randint(state, limit_ulong);
        else
            m_d = n_randint(state, limit_low_ulong);
        fmpz_pow_ui(limit_x, base, FLINT_D_BITS);
        fmpz_randm(x, state, limit_x);

        dmod_t mod;
        dmod_init(&mod, m_d);
        
        fmpz_mod_ui(x, x, mod.n);
        
        len = n_randint(state, 1000);
        
        if (!len)
            continue;

        a = _fmpz_vec_init(len);
        ans = _fmpz_vec_init(len);
        result = _fmpz_vec_init(len);
        result1 = _fmpz_vec_init(len);
        c = _dmod_vec_init(len);
        
        for (j = 0; j < len; j++)
        {
            fmpz_randtest_not_zero(a + j, state, FLINT_BITS);
            fmpz_mod_ui(a + j, a + j, mod.n);
            c[j] = fmpz_get_d(a + j);

        }
        alpha = fmpz_get_d(x);
        
        _dmod_vec_scalar_mul(c, alpha, len, mod); 
        
        _fmpz_vec_scalar_mul_fmpz(ans, a, len, x);
        
        fmpz_set_ui(mod_fmpz, mod.n);
        _fmpz_vec_scalar_mod_fmpz(result1, ans, len, mod_fmpz);
        
        for (j = 0; j < len; j++)
        {    
            fmpz_set_d(result + j, c[j]);
            
            if(fmpz_equal(result1 + j, result + j) == 0)
            {
                printf("FAIL");
                abort();
            }
        }
   
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(ans, len);
        _fmpz_vec_clear(result, len);
        _fmpz_vec_clear(result1, len);
        
        _dmod_vec_clear(c);
    
        fmpz_clear(m);
        fmpz_clear(limit);
        fmpz_clear(limit_low);
        fmpz_clear(limit_x);
        fmpz_clear(x);
        fmpz_clear(base);
        fmpz_clear(sum);
        fmpz_clear(sum1);
        fmpz_clear(mod_fmpz);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
