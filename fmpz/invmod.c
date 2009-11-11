/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

ulong z_gcdinv(ulong * inv, long a, ulong b)
{
   ulong ua = FLINT_ABS(a);
   ulong g;

   if (ua >= b) ua %= b;
   
   g = n_gcdinv(inv, ua, b);
   
   if (a < 0L) *inv = n_submod(0UL, *inv, b);

   return g;
}

int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
	fmpz c1 = *g;
	fmpz c2 = *h;
	int val;

   if (fmpz_is_zero(h))
   {
      printf("Exception : divide by zero in fmpz_invmod\n");
      abort();
   }

   if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // h is also small
		{
			ulong inv;
			if (c2 < 0L) c2 = -c2;
			if (c2 == 1L)  return 0; // special case not handled by n_invmod
			ulong gcd = z_gcdinv(&inv, c1, c2);
			if (gcd == 1UL) 
			{
				fmpz_set_si(f, inv); // check gcd is 1
				return 1;
			} else return 0;
		} else // h is large and g is small
		{
			__mpz_struct temp; // put g into a temporary mpz_t
			if (c1 < 0L) 
			{
				c1 = -c1;
				temp._mp_d = &c1;
			   temp._mp_size = -1;
			} else if (c1 == 0L) temp._mp_size = 0;
			else 
			{
				temp._mp_d = &c1;
				temp._mp_size = 1;
			}
			
			__mpz_struct * mpz_ptr = _fmpz_promote(f);
			val = mpz_invert(mpz_ptr, &temp, COEFF_TO_PTR(c2));
			_fmpz_demote_val(f); // inverse mod h may result in small value
			
			return val;
		}
	} else // g is large
	{
      if (!COEFF_IS_MPZ(c2)) // h is small
		{
			ulong inv;
			if (c2 < 0L) c2 = -c2;
			if (c2 == 1L) return 0; // special case not handled by z_gcd_invert
			// reduce g mod h first
			
			ulong r = mpz_fdiv_ui(COEFF_TO_PTR(c1), c2);
			
			ulong gcd = z_gcdinv(&inv, r, c2);
			if (gcd == 1UL) 
			{
				fmpz_set_si(f, inv); // check gcd is 1
				return 1;
			} else return 0;
		} else // both are large
		{			
			__mpz_struct * mpz_ptr = _fmpz_promote(f);
			val = mpz_invert(mpz_ptr, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2));
			_fmpz_demote_val(f); // reduction mod h may result in small value
			
			return val;
		}	
	}
}

