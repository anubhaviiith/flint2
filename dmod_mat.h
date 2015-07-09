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

#ifndef DMOD_MAT_H
#define DMOD_MAT_H

#include <math.h>
#include <gmp.h>
#include "double_extras.h"
#include "flint.h"
#include "ulong_extras.h"
#include <math.h>
#include "dmod_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif


typedef struct
{
    slong nrows;
    slong ncols;
    slong ld;
    double *rows;
    dmod_t mod;
}dmod_mat_struct;


typedef dmod_mat_struct dmod_mat_t[1];

#define dmod_mat_nrows(mat) ((mat)->nrows)
#define dmod_mat_ncols(mat) ((mat)->ncols)
#define dmod_mat_ld(mat) ((mat)->ld)

#define dmod_mat_entry(mat, i, j) ((mat)->rows[ (( (dmod_mat_ld(mat)) * (i)) +  (j)) ])
#define dmod_mat_entry_ptr(mat, i, j) ( &dmod_mat_entry(mat, i, j) )

#define DMOD_MAT_MUL_STRASSEN_CUTOFF 300

static __inline__
void _dmod_mat_set_mod(dmod_mat_t mat, double n)
{
    mat->mod.n = n;
    mat->mod.ninv = (double)1/n;
    mat->mod.nbits = FLINT_BIT_COUNT(n);
}

static __inline__
void _dmod_mat_set(dmod_mat_t mat, slong i, slong j, double val)
{
    dmod_mat_entry(mat, i, j) = val;
}

static __inline__
void _dmod_mat_print(dmod_mat_t mat)
{
    slong i, j, m, n; 
    m = dmod_mat_nrows(mat);
    n = dmod_mat_ncols(mat);
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            flint_printf("%lf ", dmod_mat_entry(mat, i, j) );
        }
        flint_printf("\n");
    }
    flint_printf("\n");
}


/*  Memory management  *******************************************************/

FLINT_DLL void _dmod_mat_init(dmod_mat_t A, slong m, slong n, dmod_t mod);

FLINT_DLL void _dmod_mat_clear(dmod_mat_t mat);



/*  Matrix-Matrix / Matrix-Vector Multiplication   *******************************************************/

FLINT_DLL void _dmod_mat_mul(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_dp(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_classical(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_strassen(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_vec(double *res, dmod_mat_t A, const double *x, slong lenx, dmod_t mod);

/*  Matrix Addition/Substraction Multiplication   *******************************************************/

FLINT_DLL void _dmod_mat_add(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_sub(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);


FLINT_DLL void _dmod_mat_randtest(dmod_mat_t A, flint_rand_t state);

FLINT_DLL void _dmod_mat_rank_update(dmod_mat_t A, const double *x, const double *y, slong lenx, slong leny, dmod_t mod);

FLINT_DLL void _dmod_mat_scalar_mul(dmod_mat_t B, const dmod_mat_t A, double c);

/* Windows *****************************************/

FLINT_DLL void _dmod_mat_window_init(dmod_mat_t window, const dmod_mat_t mat, slong r1, slong c1, slong m, slong n);

FLINT_DLL void _dmod_mat_window_clear(dmod_mat_t window);



#ifdef __cplusplus
}
#endif

#endif

