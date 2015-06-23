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

    Copyright (C) 2015 Anubhav Srivastava

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
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("mod....");
    fflush(stdout);

    dmod_t mod;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong limit_ulong, m_d, a;
        double result1, result2, c;

        limit_ulong = pow(2, FLINT_D_BITS);
        
        a = n_randint(state, limit_ulong); 
        m_d = n_randint(state, limit_ulong);
        
        dmod_init(&mod, m_d);
    
        c = (double)a;

        result1 = dmod_reduce(c, mod);
        
        result2 = n_mod2_precomp(a, mod.n, mod.ninv);

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
