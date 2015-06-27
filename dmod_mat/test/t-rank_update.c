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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"
#include "dmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, rep, rand;
    FLINT_TEST_INIT(state);
    

    flint_printf("rank-1....");
    fflush(stdout);

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        dmod_mat_t A;
        double *x, *y; 

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        
        rand = n_randint(state, 10);
            
        dmod_t mod;
        dmod_init(&mod, rand); 

        _dmod_mat_init(A, m, n, mod);
        x = _dmod_vec_init(m);  
        y = _dmod_vec_init(m);  
                
        _dmod_mat_randtest(A, state);
        _dmod_vec_randtest(x, state, m, mod);
        _dmod_vec_randtest(y, state, n, mod);
        
        _dmod_mat_rank_update(A, x, y, m, n, mod); 
        
        _dmod_mat_clear(A);
        _dmod_vec_clear(x);
        _dmod_vec_clear(y);

        break;
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
