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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_resultant(fmpz_t res, const fmpz * poly1, slong len1, 
                                               const fmpz * poly2, slong len2)
{
   slong bits1 = FLINT_ABS(_fmpz_vec_max_bits(poly1, len1)); 
   slong bits2 = FLINT_ABS(_fmpz_vec_max_bits(poly2, len2));

   if (len2 > 144 || len2*len2*len2*(bits1 + bits2) > WORD(6000000))
      _fmpz_poly_resultant_modular(res, poly1, len1, poly2, len2);
   else
      _fmpz_poly_resultant_euclidean(res, poly1, len1, poly2, len2);
}

void fmpz_poly_resultant(fmpz_t res, const fmpz_poly_t poly1,
                                                      const fmpz_poly_t poly2)
{
   slong len1 = poly1->length;
   slong len2 = poly2->length;
   
   if (len1 == 0 || len2 == 0)
     fmpz_zero(res);
   else if (len1 >= len2)
        _fmpz_poly_resultant(res, poly1->coeffs, len1, poly2->coeffs, len2);
   else
   {
        _fmpz_poly_resultant(res, poly2->coeffs, len2, poly1->coeffs, len1);  
        if ((len1 > 1) && (!(len1 & WORD(1)) & !(len2 & WORD(1))))
            fmpz_neg(res, res);
   }
}
