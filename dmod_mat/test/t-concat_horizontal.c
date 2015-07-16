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

    Copyright (C) 2015 Elena Sergeicheva

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_vec.h"
#include "dmod_mat.h"
#include "ulong_extras.h"

int main(void)
{
    dmod_mat_t A, B, C;
    dmod_mat_t window1, window2;
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("concat_horizontal....");
    fflush(stdout);


    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong c1, c2, r1, n;

        c1 = n_randint(state, 50);
        c2 = n_randint(state, 50);
        r1 = n_randint(state, 50);
        n = n_randint(state, 50) + 1;

        dmod_t mod;
        dmod_init(&mod, n); 
        
        _dmod_mat_init(A, r1, c1, mod);
        _dmod_mat_init(B, r1, c2, mod);
        _dmod_mat_init(C, r1, c1 + c2, mod);

        _dmod_mat_randtest(A, state);
        _dmod_mat_randtest(B, state);
        
        _dmod_mat_randtest(C, state);

        _dmod_mat_concat_horizontal(C, A, B);
        
        _dmod_mat_window_init(window1, C, 0, 0, r1, c1);
        _dmod_mat_window_init(window2, C, 0, c1, r1, c2);

        if (!(_dmod_mat_equal(window1, A) && _dmod_mat_equal(window2, B)))
        {
            flint_printf("A = \n");
            _dmod_mat_print(A);
            flint_printf("B = \n");
            _dmod_mat_print(B);
            flint_printf("A concat_horizontal B = \n");
            _dmod_mat_print(C);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
         
        _dmod_mat_clear(A);
        _dmod_mat_clear(B);
        _dmod_mat_clear(C);

        _dmod_mat_window_clear(window1);
        _dmod_mat_window_clear(window2);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
