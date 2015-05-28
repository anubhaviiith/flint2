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

#define DMOD_VEC_SP_EPS (1e-14)

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_addmul....");
    fflush(stdout);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        double *a, *b, *bcopy, *res;
        slong len = n_randint(state, 100);
        
        slong x = n_randint(state, 100); 
        
        if (!len)
            continue;

        a = _dmod_vec_init(len);
        b = _dmod_vec_init(len);
        bcopy = _dmod_vec_init(len);

        _dmod_vec_randtest(a, state, len, 0, 0);
        _dmod_vec_randtest(b, state, len, 0, 0);
        
        _dmod_vec_copy(b, bcopy, len);  
        
        _dmod_vec_scalar_addmul(a, b, x, len); 
 
        _dmod_vec_scalar_mul(a, x, len);

        _dmod_vec_sub(b, a, len);
        
        /*result = _dmod_vec_equal(bcopy, b, len);
        */
        slong z;
        for (z=0;z<len;++z)
        {
            if(bcopy[z]!=b[z])
            {
                printf("%lf %lf\n",bcopy[z], b[z]);
                break;
            }
        }
        if (result == 0)
        {
            flint_printf("FAIL:\n");
            abort();
        }
        
        _dmod_vec_clear(a);
        _dmod_vec_clear(b);
        _dmod_vec_clear(bcopy);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
