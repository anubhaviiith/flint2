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

int main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("dot....");
    fflush(stdout);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        double *a, *b;
        mp_limb_t res1, sum = 0;
        
        dmod_t mod;
        mp_limb_t m;

        slong len = n_randint(state, 100);
        m = n_randtest_not_zero(state);
        
        dmod_init(&mod, m);

        if (!len)
            continue;

        a = _dmod_vec_init(len);
        b = _dmod_vec_init(len);

        _dmod_vec_randtest(a, state, len, 0, 0);
        _dmod_vec_randtest(b, state, len, 0, 0);

        res1 = _dmod_vec_dot(a, b, len, mod);
        
        slong i;
        for (i = 0; i < len; i++)
        {
            sum = sum + (a[i] * b[i]);

        }
        sum = n_mod2_precomp(sum, mod.n, mod.ninv);
        
        flint_printf("%wu\n", res1);
        
        _dmod_vec_clear(a);
        _dmod_vec_clear(b);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
