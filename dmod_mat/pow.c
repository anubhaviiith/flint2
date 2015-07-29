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

    Copyright (C) 2014 Ashish Kedia

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "dmod_mat.h"

void _dmod_mat_pow(dmod_mat_t dest, const dmod_mat_t mat, ulong pow)
{
    dmod_mat_t temp1, temp2;
    
    if (mat->nrows == 0)
    {
        return;
    }
    if (pow == 0)
    {
        _dmod_mat_one(dest);
        return;
    }
    if (pow == 1)
    {
        _dmod_mat_copy(dest, mat);
        return;
    }
    
    if (pow == 2)
    {
        _dmod_mat_mul(dest, mat, mat);
        return;
    }

    _dmod_mat_init(temp1, mat->nrows, mat->ncols, mat->mod);
     
    if(pow == 3)
    {
        _dmod_mat_mul(temp1, mat, mat);
        _dmod_mat_mul(dest, temp1, mat);
        _dmod_mat_clear(temp1);
        return;
    }

    _dmod_mat_one(dest);
    _dmod_mat_init(temp2, mat->nrows, mat->ncols, mat->mod);
    _dmod_mat_copy(temp2, mat);
    
    while (pow > 0)
    {
        if (pow % 2 == 1)
        {
            _dmod_mat_mul(temp1, dest, temp2);
            _dmod_mat_swap(temp1, dest);
        }
        if (pow > 1)
        {
            _dmod_mat_mul(temp1, temp2, temp2);
            _dmod_mat_swap(temp1, temp2);
        }
        pow /= 2;
    }
    
    _dmod_mat_clear(temp1);
    _dmod_mat_clear(temp2);
}
void dmod_mat_pow(dmod_mat_t dest, const dmod_mat_t mat, ulong pow)
{
    dmod_mat_t temp;
    
    if (mat->ncols == 0 || mat->nrows == 0)
        return;

    if (mat == dest)
    {
        _dmod_mat_init(temp, mat->nrows, mat->ncols, mat->mod);
        _dmod_mat_copy(temp, mat);
        _dmod_mat_pow(dest, temp, pow);
        _dmod_mat_clear(temp);
    }
    else
    {
        _dmod_mat_pow(dest, mat, pow);
    }
}
