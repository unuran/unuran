/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      gamma.c                                                      *
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

#include <unur_methods.h>
#include <unur_distribution.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
/* Prototypes for special generators                                         */

inline static double gll(double a, UNUR_URNG_TYPE urng);
/* Rejection with log-logistic envelopes                                     */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define GEN     gen->data.cstd
#define uniform()  (_unur_call_urng_prt(urng))

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  get special sampling routine for distribution                          **/
/**                                                                         **/
/*****************************************************************************/

_UNUR_SAMPLING_ROUTINE_CONT *
_unur_stdgen_gamma_get_routine(unsigned variant)
     /*----------------------------------------------------------------------*/
     /* get pointer to sampling routine                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   variant ... variant of special generator                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to sampling routine                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  switch (variant) {
  case 0: /* Default */
    return unur_stdgen_sample_gamma_gll; /* Rejection with log-logistic envelopes */
  case UNUR_STDGEN_INVERSION:
  default:
    return NULL;
  }

} /* end of _unur_stdgen_gamma_get_routine() */

/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO

const char *
_unur_stdgen_gamma_routinename(void *routine)
     /*----------------------------------------------------------------------*/
     /* get name of sampling routine                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   routine ... pointer to sampling routine                            */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to name of sampling routine                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
#define routinename(rn) if (routine==(rn)) return #rn

  routinename( unur_stdgen_sample_gamma_gll );

  return NULL;

} /* end of _unur_stdgen_gamma_routinename() */

#endif

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Wrapper for special generators                                         **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_gamma_gll( struct unur_gen *gen )
     /* Rejection with log-logistic envelopes                                */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (gll( GEN.pdf_param[0],gen->urng ) * GEN.pdf_param[1] + GEN.pdf_param[2]);

} /* end of unur_stdgen_sample_gamma_gll() */


/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators (WinRand)                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

inline double
gll(double a, UNUR_URNG_TYPE urng)
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
{
 static double aa,bb,cc,a_in = -1.0;
 double u1,u2,v,r,z,gl;

 if (a != a_in) {
   a_in = a;
   aa = (a > 1.0) ? sqrt(a + a - 1.0) : a;
   bb = a - 1.386294361;
   cc = a + aa;
 }
 while (1) {
   u1 = uniform();
   u2 = uniform();
   v = log(u1 / (1.0 - u1)) / aa;
   gl = a * exp(v);
   r = bb + cc * v - gl;
   z = u1 * u1 * u2;
   if (r + 2.504077397 >= 4.5 * z) break;
   if (r >= log(z)) break;
 }
 return gl;
} /* end of gll() */

/*---------------------------------------------------------------------------*/
