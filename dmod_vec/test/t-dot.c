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

#define EPSILON 0.0000000000001

int main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("dot....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        double *a, *b;
        double m, result1 = 0, result2 = 0; 
        slong len = n_randint(state, 1000);
        
        m = n_randint(state, n_powmod_precomp(2,(FLINT_D_BITS - FLINT_FLOG2(len))/2 , 1, 1) - 1);
 
        dmod_t mod;
        dmod_init(&mod, m);

        if (!len)
            continue;

        a = _dmod_vec_init(len);
        b = _dmod_vec_init(len);

        _dmod_vec_randtest(a, state, len, 0, 0);
        _dmod_vec_randtest(b, state, len, 0, 0);

        /*Code to be tested */
        result1 = _dmod_vec_dot(a, b, len, mod); /* returns double */
        
        /*Test */

        for (j = 0; j < len; ++j)
        {
            a[j] = dmod_precomp(a[j], mod.n, mod.ninv);
            b[j] = dmod_precomp(b[j], mod.n, mod.ninv);
            if ( j != 0 && j % (n_powmod_precomp(2,(FLINT_D_BITS - (2 * mod.b)), mod.n, mod.ninv)) == 0)
            {
                result2 = dmod_precomp(result2, mod.n, mod.ninv);
            }
            result2 += a[j]*b[j];
        }
        
        if ( fabs(result1 - result2) > EPSILON)
        {
            flint_printf("FAIL\n");
            abort();
        }

        _dmod_vec_clear(a);
        _dmod_vec_clear(b);
        
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
