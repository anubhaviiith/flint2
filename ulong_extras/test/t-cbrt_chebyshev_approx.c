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

    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);
   
    flint_printf("cbrt_chebyshev_approx....");
    fflush(stdout);

    /* random n */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);
      
        n = n_randtest(state);
        val = n_cbrt_chebyshev_approx(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);
      
        result = (val == ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans); 
            abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }
    /* type n^3 + k */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans, bits;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);
      
        bits = n_randint(state, FLINT_BITS/3 + 1);
        n = n_randtest_bits(state, bits);
        n = n*n*n;
        n += (n_randint(state, 100) - 50);
        val = n_cbrt_chebyshev_approx(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);
      
        result = (val == ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans); 
            abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }

    /* type n^3 + 1 */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans, bits;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);
        
        bits = n_randint(state, FLINT_BITS/3 + 1);
        n = n_randtest_bits(state, bits);
        n = n*n*n;
        n += 1;
        val = n_cbrt_chebyshev_approx(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);
      
        result = (val == ans);
        if (!result)
        { 
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans); 
            abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }
   
    /* type n^3 - 1 */

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n, val, ans, bits;
        mpz_t mpz_n, mpz_val;

        mpz_init(mpz_n);
        mpz_init(mpz_val);
      
        bits = n_randint(state, FLINT_BITS/3 + 1);
        n = n_randtest_bits(state, bits);
        n = n*n*n;
        n -= 1;
        val = n_cbrt_chebyshev_approx(n);

        flint_mpz_set_ui(mpz_n, n);
        mpz_root(mpz_val, mpz_n, 3);
        ans = flint_mpz_get_ui(mpz_val);
      
        result = (val == ans);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, val = %wd, ans = %wu\n", n, val, ans); 
            abort();
        }
        mpz_clear(mpz_n);
        mpz_clear(mpz_val);
    }

   FLINT_TEST_CLEANUP(state);
   flint_printf("PASS\n");
   return 0;
}
