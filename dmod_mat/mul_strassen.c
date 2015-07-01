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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_vec.h"
#include "dmod_mat.h"


void _dmod_mat_mul_strassen(dmod_mat_t C, const dmod_mat_t A, const dmod_mat_t B)
{
    #if HAVE_BLAS
    slong a, b, c;
    slong anr, anc, bnr, bnc;

    dmod_mat_t A11, A12, A21, A22;
    dmod_mat_t B11, B12, B21, B22;
    dmod_mat_t C11, C12, C21, C22; 
    dmod_mat_t X1, X2;

    a = A->nrows;
    b = A->ncols;
    c = B->ncols;

    if (a <= 4 || b <= 4 || c <= 4)
    {
        _dmod_mat_mul(C, A, B);
        return;
    }

    anr = a / 2;
    anc = b / 2;
    bnr = anc;
    bnc = c / 2;
    
    _dmod_mat_window_init(A11, A, 0, 0, anr, anc);
    _dmod_mat_window_init(A12, A, 0, anc, anr, anc);
    _dmod_mat_window_init(A21, A, anr, 0, anr, anc);
    _dmod_mat_window_init(A22, A, anr, anc, anr, anc);

    _dmod_mat_window_init(B11, B, 0, 0, bnr, bnc);
    _dmod_mat_window_init(B12, B, 0, bnc, bnr, bnc);
    _dmod_mat_window_init(B21, B, bnr, 0, bnr, bnc);
    _dmod_mat_window_init(B22, B, bnr, bnc, bnr, bnc);

    _dmod_mat_window_init(C11, C, 0, 0, anr, bnc);
    _dmod_mat_window_init(C12, C, 0, bnc, anr, bnc);
    _dmod_mat_window_init(C21, C, anr, 0, anr, bnc);
    _dmod_mat_window_init(C22, C, anr, bnc, anr, bnc);

    _dmod_mat_init(X1, anr, FLINT_MAX(bnc, anc), A->mod);
    _dmod_mat_init(X2, anc, bnc, A->mod);

    X1->ncols = anc;
    X1->ld = anc;

    ulong zero = 0;

    _dmod_mat_sub(X1, A11, A21); 
    _dmod_mat_sub(X2, B22, B12);
    _dmod_mat_mul(C21, X1, X2);

    _dmod_mat_add(X1, A21, A22);
    _dmod_mat_sub(X2, B12, B11);
    _dmod_mat_mul(C22, X1, X2);

    _dmod_mat_sub(X1, X1, A11);
    _dmod_mat_sub(X2, B22, X2);
    _dmod_mat_mul(C12, X1, X2);

    _dmod_mat_sub(X1, A12, X1);
    _dmod_mat_mul(C11, X1, B22);
    

    X1->ncols = bnc;
    X1->ld = bnc;

    _dmod_mat_mul(X1, A11, B11);

    _dmod_mat_add(C12, X1, C12);
    _dmod_mat_add(C21, C12, C21);
    _dmod_mat_add(C12, C12, C22);
    _dmod_mat_add(C22, C21, C22);
    _dmod_mat_add(C12, C12, C11);
    _dmod_mat_sub(X2, X2, B21);
    _dmod_mat_mul(C11, A22, X2);
    _dmod_mat_sub(C21, C21, C11);
    _dmod_mat_mul(C11, A12, B21);
    _dmod_mat_add(C11, X1, C11);


    _dmod_mat_clear(X1);
    _dmod_mat_clear(X2);

    _dmod_mat_window_clear(A11);
    _dmod_mat_window_clear(A12);
    _dmod_mat_window_clear(A21);
    _dmod_mat_window_clear(A22);

    _dmod_mat_window_clear(B11);
    _dmod_mat_window_clear(B12);
    _dmod_mat_window_clear(B21);
    _dmod_mat_window_clear(B22);

    _dmod_mat_window_clear(C11);
    _dmod_mat_window_clear(C12);
    _dmod_mat_window_clear(C21);
    _dmod_mat_window_clear(C22);


    if (c > 2*bnc) 
    {
        dmod_mat_t Bc, Cc;
        _dmod_mat_window_init(Bc, B, 0, 2*bnc, b, c - 2*bnc);
        _dmod_mat_window_init(Cc, C, 0, 2*bnc, a, c - 2*bnc);
        _dmod_mat_mul(Cc, A, Bc);
        _dmod_mat_window_clear(Bc);
        _dmod_mat_window_clear(Cc);
    }

    if (a > 2*anr) 
    {
        dmod_mat_t Ar, Cr;
        _dmod_mat_window_init(Ar, A, 2*anr, 0, a - 2*anr, b);
        _dmod_mat_window_init(Cr, C, 2*anr, 0, a - 2*anr, c);
        _dmod_mat_mul(Cr, Ar, B);
        _dmod_mat_window_clear(Ar);
        _dmod_mat_window_clear(Cr);
    }

    if (b > 2*anc) 
    {
        dmod_mat_t Ac, Br, Cb, temp;

        _dmod_mat_window_init(Ac, A, 0, 2*anc, 2*anr, b - 2*anc);
        _dmod_mat_window_init(Br, B, 2*bnr, 0, b - 2*bnr, 2*bnc);
        _dmod_mat_window_init(Cb, C, 0, 0, 2*anr, 2*bnc);
        _dmod_mat_init(temp, 2*anr, 2*bnc, A->mod);

        _dmod_mat_mul(temp, Ac, Br); 
        _dmod_mat_add(Cb, Cb, temp);

        _dmod_mat_clear(temp);

        _dmod_mat_window_clear(Ac);
        _dmod_mat_window_clear(Br);
        _dmod_mat_window_clear(Cb);

   }
   #endif
}
