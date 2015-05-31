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

#ifndef DMOD_VEC_H
#define DMOD_VEC_H

#include <math.h>
#include <gmp.h>
#include "double_extras.h"
#include "flint.h"
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
   mp_limb_t n;
   mp_limb_t ninv;
} dmod_t;


static __inline__
void dmod_init(dmod_t * mod, mp_limb_t n)
{
   mod->n = n;
   mod->ninv = n_precompute_inverse(n);
}

/*  Memory management  *******************************************************/

FLINT_DLL double * _dmod_vec_init(slong len);

FLINT_DLL void _dmod_vec_clear(double * vec);

/*  Randomisation  ***********************************************************/

FLINT_DLL void _dmod_vec_randtest(double * f, flint_rand_t state, slong len, slong minexp, slong maxexp);

/*  Dot product  **************************************/

FLINT_DLL mp_limb_t _dmod_vec_dot(const double * vec1, const double * vec2, slong len2, dmod_t mod);

/* Substraction **************************************/

FLINT_DLL void  _dmod_vec_sub(double * vec1, const double * vec2, slong len2);

/* Is equal  **************************************/

FLINT_DLL int  _dmod_vec_equal(const double * vec1, const double * vec2, slong len2);

/* Scalar mul and scalar addmul **************************************/

FLINT_DLL void _dmod_vec_scalar_mul(double * vec1, const double alpha, slong len2);

FLINT_DLL void _dmod_vec_scalar_addmul(const double * vec1, double * vec2, const double alpha, slong len2);

/* Copy  **************************************/

FLINT_DLL void _dmod_vec_copy(const double * vec1, double * vec2, slong len2);

#ifdef __cplusplus
}
#endif

#endif
