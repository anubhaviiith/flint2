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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

double d_invmod1(double x, double y)
{
    double v1 = 0;
    double v2 = 1;
    double t2;
    double u3, v3;
    ulong quot, rem;
    
    u3 = y, v3 = x;

    if (v3 > u3)
    {
        rem = u3;
        u3 = v3;
        t2 = v2;
        v2 = v1;
        v1 = t2;
        v3 = rem;
    }

    while (v3)
    {
        quot = u3 / v3;
        rem = u3 - v3 * quot;
        u3 = v3;
        t2 = v2;
        v2 = v1 - quot * v2;
        v1 = t2;
        v3 = rem;
    }

    if (v1 < 0.0)
        v1 += y;
     
    return v1;
}
