/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      gamma_gen.c                                                  *
 *                                                                           *
 *   Special generators for Gamma distribution                               *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_defs.h>

#include <unur_methods.h>
#include <unur_distribution.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static void gamma_gll_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.cstd        /* data for parameter object         */
#define GEN       gen->data.cstd        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

#define alpha (DISTR.params[0])
#define beta  (DISTR.params[1])
#define gamma (DISTR.params[2])

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_gamma_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for gamma distribution                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0.);
  COOKIE_CHECK(par,CK_CSTD_PAR,0.);

  switch (par->variant) {
  case 0:  /* Rejection with log-logistic envelopes */  /* DEFAULT */
    _unur_cstd_set_sampling_routine( par,gen,unur_stdgen_sample_gamma_gll );
    gamma_gll_init( gen );
    return 1;
  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default:
    return 0;
  }

} /* end of _unur_stdgen_gamma_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Gamma Distribution: Rejection from log-logistic envelopes                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               Gamma distribution with parameter a.                        *
 *                                                                           *
 * REFERENCE:  - R.C.H. Cheng (1977): The Generation of Gamma Variables      *
 *               with Non-Integral Shape Parameter,                          *
 *               Appl. Statist. 26(1), 71-75                                 *
 *                                                                           *
 * Implemented by Ralf Kremer                                                *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

#define aa  GEN.gen_param[0]
#define bb  GEN.gen_param[1]
#define cc  GEN.gen_param[2]

inline static void
gamma_gll_init( struct unur_gen *gen )
{
  if (GEN.gen_param == NULL) {
    GEN.gen_param = _unur_malloc(3 * sizeof(double));
    GEN.n_gen_param = 3;
  }

  /* -X- setup code -X- */
  aa = (alpha > 1.0) ? sqrt(alpha + alpha - 1.0) : alpha;
  bb = alpha - 1.386294361;
  cc = alpha + aa;
  /* -X- end of setup code -X- */

} /* end of gamma_gll_init() */

double 
unur_stdgen_sample_gamma_gll( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X;
  double u1,u2,v,r,z;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  while (1) {
    u1 = uniform();
    u2 = uniform();
    v = log(u1 / (1.0 - u1)) / aa;
    X = alpha * exp(v);
    r = bb + cc * v - X;
    z = u1 * u1 * u2;
    if (r + 2.504077397 >= 4.5 * z) break;
    if (r >= log(z)) break;
  }
  /* -X- end of generator code -X- */

  return (X * beta + gamma);

} /* end of unur_stdgen_sample_gamma_gll() */

#undef aa
#undef bb
#undef cc

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
#undef gamma
/*---------------------------------------------------------------------------*/
