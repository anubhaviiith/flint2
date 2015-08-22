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

#define DMOD_MAT_MUL_STRASSEN_CUTOFF 5000
#define DMOD_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define DMOD_MAT_SOLVE_TRI_COLS_CUTOFF 64
#define DMOD_MAT_LU_RECURSIVE_CUTOFF 4

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

/* Comparisons and assignments/swap  ***********************************************/

FLINT_DLL void _dmod_mat_copy(dmod_mat_t B, const dmod_mat_t A);

FLINT_DLL int _dmod_mat_equal(const dmod_mat_t mat1, const dmod_mat_t mat2);

FLINT_DLL void _dmod_mat_one(dmod_mat_t A);

FLINT_DLL void _dmod_mat_swap(dmod_mat_t mat1, dmod_mat_t mat2);


/*  Matrix-Matrix / Matrix-Vector Multiplication   *******************************************************/

FLINT_DLL void _dmod_mat_mul(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_dp(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_classical(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_strassen(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_mul_vec(double *res, dmod_mat_t A, const double *x, slong lenx, dmod_t mod);

/*  Arithmetic Operations   *******************************************************/

FLINT_DLL void _dmod_mat_add(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_sub(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_scalar_mul(dmod_mat_t B, const dmod_mat_t A, double c);

FLINT_DLL double dmod_mat_det(const dmod_mat_t A);

FLINT_DLL double _dmod_mat_trace(const dmod_mat_t mat);

FLINT_DLL int _dmod_mat_inv(dmod_mat_t B, const dmod_mat_t A);

FLINT_DLL void dmod_mat_pow(dmod_mat_t dest, const dmod_mat_t mat, ulong pow);

FLINT_DLL slong _dmod_mat_nullspace(dmod_mat_t X, const dmod_mat_t A);

FLINT_DLL slong _dmod_mat_rref(dmod_mat_t A, slong * pivots_nonpivots, slong * P);

FLINT_DLL void _dmod_mat_transpose(dmod_mat_t B, const dmod_mat_t A);

/*Random ******************************/

FLINT_DLL void _dmod_mat_randtest(dmod_mat_t A, flint_rand_t state);

/*Concatenation ***********************/

FLINT_DLL void _dmod_mat_concat_vertical(dmod_mat_t res, const dmod_mat_t mat1, const dmod_mat_t mat2);

FLINT_DLL void _dmod_mat_concat_horizontal(dmod_mat_t res, const dmod_mat_t mat1, const dmod_mat_t mat2);


/*LU decomposition and rank **********************/

FLINT_DLL slong _dmod_mat_lu(slong * P, dmod_mat_t A, int rank_check);

FLINT_DLL slong _dmod_mat_lu_classical(slong * P, dmod_mat_t A, int rank_check);

FLINT_DLL slong _dmod_mat_lu_recursive(slong * P, dmod_mat_t A_d, int rank_check);

FLINT_DLL slong _dmod_mat_rank(const dmod_mat_t A);


/*Reduction *****************/

FLINT_DLL void _dmod_mat_reduce(dmod_mat_t C);


/*Solve ********************/

FLINT_DLL int _dmod_mat_solve(dmod_mat_t X, const dmod_mat_t A, const dmod_mat_t B);

FLINT_DLL void _dmod_mat_solve_triu(dmod_mat_t X, const dmod_mat_t L, const dmod_mat_t B, int unit);

FLINT_DLL void _dmod_mat_solve_triu_classical(dmod_mat_t X, const dmod_mat_t U, const dmod_mat_t B, int unit);

FLINT_DLL void _dmod_mat_solve_triu_recursive(dmod_mat_t X, const dmod_mat_t U, const dmod_mat_t B, int unit);

FLINT_DLL void _dmod_mat_solve_tril(dmod_mat_t X, const dmod_mat_t L, const dmod_mat_t B, int unit);

FLINT_DLL void _dmod_mat_solve_tril_classical(dmod_mat_t X, const dmod_mat_t L, const dmod_mat_t B, int unit);

FLINT_DLL void _dmod_mat_solve_tril_recursive(dmod_mat_t X, const dmod_mat_t L, const dmod_mat_t B, int unit);

/*Window **********************/

FLINT_DLL void _dmod_mat_window_init(dmod_mat_t window, const dmod_mat_t A, slong r1, slong c1, slong m, slong n);

FLINT_DLL void _dmod_mat_window_clear(dmod_mat_t mat);

#ifdef __cplusplus
}
#endif

#endif

