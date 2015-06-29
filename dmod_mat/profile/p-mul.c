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
#include "nmod_mat.h"
#include "dmod_mat.h"
#include "dmod_vec.h"
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
    ulong i, j, dim = params->dim;
    int algorithm = params->algorithm;

    nmod_mat_t A, B, C;
    dmod_mat_t A_d, B_d, C_d;

    nmod_mat_init(A, dim, dim, n);
    nmod_mat_init(B, dim, dim, n);
    nmod_mat_init(C, dim, dim, n);
    
    dmod_t mod;
    dmod_init(&mod, n);

    _dmod_mat_init(A_d, dim, dim, mod);
    _dmod_mat_init(B_d, dim, dim, mod);
    _dmod_mat_init(C_d, dim, dim, mod);
    
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            _dmod_mat_set(A_d, i, j, (double)A->rows[i][j]);
            _dmod_mat_set(B_d, i, j, (double)B->rows[i][j]);
            _dmod_mat_set(C_d, i, j, (double)C->rows[i][j]);
        }
    }
   
    prof_start();

    if (algorithm == 1)
        for (i = 0; i < count; i++)
            _dmod_mat_mul(C_d, A_d, B_d);
    else if (algorithm == 2)
        for (i = 0; i < count; i++)
            nmod_mat_mul(C, A, B);

    prof_stop();

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);

    _dmod_mat_clear(A_d);
    _dmod_mat_clear(B_d);
    _dmod_mat_clear(C_d);

}

int main(void)
{
    double min_blas, min_noblas, max;
    mat_mul_t params;
    slong dim;

    flint_printf("dmod_mat_mul:\n");
    
    params.modulus = 4000000;

    for (dim = 2; dim <= 1200; dim = (slong) ((double) dim * 1.1) + 1)
    {
        params.dim = dim;

        params.algorithm = 1;
        prof_repeat(&min_blas, &max, sample, &params);

        params.algorithm = 2;
        prof_repeat(&min_noblas, &max, sample, &params);

        flint_printf("dim = %wd, BLAS %.2f us NON-BLAS %.2f us\n", dim, min_blas, min_noblas);
    }

    return 0;
}
