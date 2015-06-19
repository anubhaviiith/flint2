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
    #if HAVE_BLAS
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("add....");
    fflush(stdout);

    dmod_t mod;
    nmod_t modn;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mp_ptr a, b, result2;
        double *c, *d, *result1;

        mp_limb_t limit_ulong, m_d, len;
       
        limit_ulong = pow(2, FLINT_D_BITS - 1);
        m_d = n_randint(state, limit_ulong);
        
        dmod_init(&mod, m_d);
        nmod_init(&modn, m_d);
        
        len = n_randint(state, 1000);
        
        if (!len)
            continue;

        a = _nmod_vec_init(len);
        b = _nmod_vec_init(len);
        result2 = _nmod_vec_init(len);

        c = _dmod_vec_init(len);
        d = _dmod_vec_init(len);
        result1 = _dmod_vec_init(len);
        
        _nmod_vec_randtest(a, state, len, modn);
        _nmod_vec_randtest(b, state, len, modn);
        
        for (j = 0; j < len; j++)
        {
            c[j] = (double)a[j];
            d[j] = (double)b[j];
        }
        
        _nmod_vec_add(result2, a, b, len, modn);

        _dmod_vec_add(result1, c, d, len, mod); 
        
       
        for (j = 0; j < len; j++)
        {    
            if(result1[j] != (double)result2[j])
            {
                printf("FAIL");
                abort();
            }
        }
   
        _dmod_vec_clear(c);
        _dmod_vec_clear(d);
        _dmod_vec_clear(result1);

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(result2);
   }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
    #endif
}
