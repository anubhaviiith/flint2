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
    slong a, b, c;
    slong anr, anc, bnr, bnc;

    dmod_mat_t X1, X2;

    a = A->nrows;
    b = A->ncols;
    c = B->ncols;

    anr = a / 2;
    anc = b / 2;
    bnr = anc;
    bnc = c / 2;

    _dmod_mat_init(X1, anr, FLINT_MAX(bnc, anc), A->mod);
    _dmod_mat_init(X2, anc, bnc, A->mod);

    X1->ncols = anc;
    X1->ld = anc;

    ulong zero = 0;

    _dmod_mat_sub_window(X1, A, A, zero, zero, zero, zero, anr, zero, anr, anc); 
    _dmod_mat_sub_window(X2, B, B, zero, zero, bnr, bnc, zero, bnc, anc, bnc);
    _dmod_mat_mul_window(C, X1, X2, anr, zero, zero, zero, zero ,zero, anr, bnc, anc);
    
    _dmod_mat_add_window(X1, A, A, zero, zero, anr, zero, anr, anc, anr, anc);
    _dmod_mat_sub_window(X2, B, B, zero, zero, zero, bnc, zero, zero, anc, bnc);
    _dmod_mat_mul_window(C, X1, X2, anr, bnc, zero, zero, zero, zero, anr, bnc, anc);

    _dmod_mat_sub_window(X1, X1, A, zero, zero, zero, zero, zero, zero, anr, anc);
    _dmod_mat_sub_window(X2, B, X2, zero, zero, bnr, bnc, zero, zero, anc, bnc);
    _dmod_mat_mul_window(C, X1, X2, zero, bnc, zero, zero, zero, zero, anr, bnc, anc);

    _dmod_mat_sub_window(X1, A, X1, zero, zero, zero, anc, zero, zero, anr, anc);
    _dmod_mat_mul_window(C, X1, B, zero, zero, zero, zero, bnr, bnc, anr, bnc, anc);
    

    X1->ncols = bnc;
    X1->ld = bnc;

    _dmod_mat_mul_window(X1, A, B, zero, zero, zero, zero, zero, zero, anr, bnc, anc);

    _dmod_mat_add_window(C, X1, C, zero, bnc, zero, zero, zero, bnc, anr, bnc);
    _dmod_mat_add_window(C, C, C, anr, zero, zero, bnc, anr, zero, anr, bnc);
    _dmod_mat_add_window(C, C, C, zero, bnc, zero, bnc, anr, bnc, anr, bnc);
    _dmod_mat_add_window(C, C, C, anr, bnc, anr, zero, anr, bnc, anr, bnc);
    _dmod_mat_add_window(C, C, C, zero, bnc, zero, bnc, zero, zero, anr, bnc);
    _dmod_mat_sub_window(X2, X2, B, zero, zero, zero, zero, bnr, zero, anc, bnc);
    _dmod_mat_mul_window(C, A, X2, zero, zero, anr, anc, zero, zero, anr, bnc, anc);
    _dmod_mat_sub_window(C, C, C, anr, zero, anr, zero, zero, zero, anr, bnc); 
    _dmod_mat_mul_window(C, A, B, zero, zero, zero, anc, bnr, zero, anr, bnc, bnr);
    _dmod_mat_add_window(C, X1, C, zero, zero, zero, zero, zero, zero, anr, bnc);


    _dmod_mat_clear(X1);
    _dmod_mat_clear(X2);

    if (c > 2*bnc) 
    {
       _dmod_mat_mul_window(C, A, B, zero, 2*bnc, zero, zero, zero, 2*bnc, a, (c - 2*bnc), b);
    }

    if (a > 2*anr) 
    {
        _dmod_mat_mul_window(C, A, B, 2*anr, zero, 2*anr, zero, zero, zero, (a - 2*anr), c, b);
    }

    if (b > 2*anc) 
    {
        dmod_mat_t temp;

        _dmod_mat_init(temp, 2*anr, 2*bnc, A->mod);
        _dmod_mat_mul_window(temp, A, B, zero, zero, zero, 2*anc, 2*bnr, zero, 2*anr, 2*bnc, (b - 2*bnr)); 
        _dmod_mat_add_window(C, C, temp, zero, zero, zero, zero, zero, zero, 2*anr, 2*bnc);
        
        _dmod_mat_clear(temp);
    }
}
