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
#include "profiler.h"
#include "flint.h"
#include "dmod_mat.h"
#include "nmod_mat.h"
#include "dmod_vec.h"
#include <omp.h>
#include "ulong_extras.h"

typedef struct
{
    ulong dim;
    mp_limb_t modulus;
    int algorithm;
} mat_mul_t;

void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    mp_limb_t n = params->modulus;
    ulong i, j, dim = params->dim, rand, q, w;
    int algorithm = params->algorithm;

    nmod_mat_t A, LU;
    dmod_mat_t A_d, LU_d;

    slong * P;
    FLINT_TEST_INIT(state);
    slong bits = n_randint(state, 27);
    if (bits < 2)
        bits = 2;
    rand = n_randprime(state, bits, 0);
    dmod_t mod;
    dmod_init(&mod, rand);  
    
    ulong r = n_randint(state, dim);

    nmod_mat_init(A, dim, dim, rand);
    nmod_mat_randrank(A, state, r);

    nmod_mat_init_set(LU, A);

    _dmod_mat_init(A_d, dim, dim, mod);
    _dmod_mat_init(LU_d, dim, dim, mod);

    for (q = 0; q < dim; q++)
    {
        for (w = 0; w < dim; w++)
        {
            dmod_mat_entry(A_d, q, w) = (double)A->rows[q][w];
        }
    }
    _dmod_mat_copy(LU_d, A_d);

    prof_start();
 
    if (algorithm == 1)
    {  
        for (i = 0; i < count; i++)
        {
            P = flint_malloc(sizeof(slong) * dim);
            _dmod_mat_lu_recursive(P, LU_d, 0); 
            flint_free(P);
        }
    }
    else if (algorithm == 2)
    {
        for (i = 0; i < count; i++)
        {
            P = flint_malloc(sizeof(slong) * dim);
            _dmod_mat_lu_classical(P, LU_d, 0);
            flint_free(P);
        }
    }
    else if (algorithm == 3)
    {
        for (i = 0; i < count; i++)
        {
            P = flint_malloc(sizeof(slong) * dim);
            nmod_mat_lu_recursive(P, LU, 0);
            flint_free(P);
        }
    }
    prof_stop();

    nmod_mat_clear(A);
    nmod_mat_clear(LU);

    _dmod_mat_clear(A_d);
    _dmod_mat_clear(LU_d);

}

int main(void)
{
    double nmodmul, dmodmul, dmodmul1, nmodstrassen, dmodclassical, dmodstrassen, max;
    mat_mul_t params;
    slong dim;

    flint_printf("dmod_mat_lu:\n");

    params.modulus = 107;

    for (dim = 2; dim <= 5000; dim = (slong) ((double) dim * 1.5) + 1)
    {
        params.dim = dim;

        params.algorithm = 1;
        prof_repeat(&nmodmul, &max, sample, &params);

        params.algorithm = 2;
        prof_repeat(&dmodmul, &max, sample, &params);
    
        params.algorithm = 3;
        prof_repeat(&dmodmul1, &max, sample, &params);


        flint_printf("%wd %.2f %.2f %0.2f\n", dim, nmodmul, dmodmul, dmodmul1);
    }


    return 0;
}
