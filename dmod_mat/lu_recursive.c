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


    Loosely based on the recursive PLS implementation in M4RI,
    Copyright (C) 2008 Clement Pernet.

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include "dmod_mat.h"

static void _apply_permutation(slong * AP, dmod_mat_t A, slong * P, slong n, slong offset)
{
    if (n != 0)
    {
        slong * APtmp;
        slong i,j;

        dmod_mat_t Atemp;
        _dmod_mat_init(Atemp, n, A->ncols, A->mod);

        APtmp = flint_malloc(sizeof(slong) * n);

        for (i = 0; i < n; i++) 
            _dmod_vec_copy(dmod_mat_entry_ptr(A, P[i] + offset, 0), dmod_mat_entry_ptr(Atemp, i, 0), A->ncols);
        
        for (i = 0; i < n; i++) 
            _dmod_vec_copy(dmod_mat_entry_ptr(Atemp, i, 0), dmod_mat_entry_ptr(A, i + offset, 0), A->ncols);

        for (i = 0; i < n; i++) APtmp[i] = AP[P[i] + offset];
        for (i = 0; i < n; i++) AP[i + offset] = APtmp[i];

        _dmod_mat_clear(Atemp); 
        flint_free(APtmp);
    }
}

slong _dmod_mat_lu_recursive(slong * P, dmod_mat_t A_d, int rank_check)
{
    slong i, j, m, n, r1, r2, n1;

    m = A_d->nrows;
    n = A_d->ncols;

    dmod_mat_t A00, A01, A10, A11, A0, A1;

    slong *P1, *P2;
  

    if (m < DMOD_MAT_LU_RECURSIVE_CUTOFF || n < DMOD_MAT_LU_RECURSIVE_CUTOFF)
    {
        r1 = _dmod_mat_lu_classical(P, A_d, rank_check);
        return r1;
    }
    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    P1 = flint_malloc(sizeof(slong) * m);
    P2 = flint_malloc(sizeof(slong) * m);

    _dmod_mat_window_init(A0, A_d, 0, 0, m, n1);
    _dmod_mat_window_init(A1, A_d, 0, n1, m, n - n1);

    r1 = _dmod_mat_lu(P1, A0, rank_check);

    if (rank_check && (r1 != n1))
    {
        flint_free(P1);
        _dmod_mat_window_clear(A0);
        _dmod_mat_window_clear(A1);
        return 0;
    }
   
    if (r1 != 0)
    {
        _apply_permutation(P, A1, P1, m, 0);
    }

    _dmod_mat_window_init(A00, A_d, 0, 0, r1, r1);
    _dmod_mat_window_init(A10, A_d, r1, 0, m - r1, r1);
    _dmod_mat_window_init(A01, A_d, 0, n1, r1, n - n1);
    _dmod_mat_window_init(A11, A_d, r1, n1, m - r1, n - n1);
    
       
    if (r1 != 0)
    {
        dmod_mat_t tmp;
        _dmod_mat_init(tmp, m - r1, n - n1, A_d->mod);
     
        _dmod_mat_solve_tril(A01, A00, A01, 1);
        _dmod_mat_mul(tmp, A10, A01);
        _dmod_mat_sub(A11, A11, tmp);
        
        _dmod_mat_clear(tmp); 
    }
    
    r2 = _dmod_mat_lu(P1, A11, rank_check);

    for (i = 0; i < m - r1; i++)
        P2[P1[i]] = i;
    
    _apply_permutation(P2, A11, P2, m - r1, 0);
    
    if (rank_check && (r1 + r2 < FLINT_MIN(m, n)))
    {
        r1 = r2 = 0;
    }
    else
    {
        _apply_permutation(P, A_d, P1, m - r1, r1);
        
        if (r1 != n1)
        {
            for (i = 0; i < m - r1; i++)
            {
                double *row = dmod_mat_entry_ptr(A_d, r1 + i, 0);
                for (j = 0; j < FLINT_MIN(i, r2); j++)
                {
                    row[r1 + j] = row[n1 + j];
                    row[n1 + j] = 0;
                }
            }
        }

        
    }
   
    flint_free(P1);
    flint_free(P2);
    _dmod_mat_window_clear(A00);
    _dmod_mat_window_clear(A01);
    _dmod_mat_window_clear(A10);
    _dmod_mat_window_clear(A11);
    _dmod_mat_window_clear(A0);
    _dmod_mat_window_clear(A1);
    
    return r1 + r2;
}
