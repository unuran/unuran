/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      powerexponential_gen.c                                       *
 *                                                                           *
 *   Special generators for Power-exponential distribution                   *
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

#include <source_unuran.h>
#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static void powerexponential_epd_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.cstd        /* data for parameter object         */
#define GEN       gen->data.cstd        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

/*---------------------------------------------------------------------------*/
#define delta  (DISTR.params[0])        /* shape */
#define theta  (DISTR.params[1])        /* location */
#define phi    (DISTR.params[2])        /* scale */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_powerexponential_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for power-exponential distribution      */
     /* if gen == NULL then only check existance of variant.                 */
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

  case 0:  /* DEFAULT */
  case 1:  /* Transformed density rejection */
    _unur_cstd_set_sampling_routine( par,gen,unur_stdgen_sample_powerexponential_epd );
    powerexponential_epd_init( gen );
    return 1;

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(par->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }

} /* end of _unur_stdgen_powerexponential_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Powerexponential Distribution: Transformed density rejection              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Power-exponential          *
 *               distribution with parameter delta <= 2 by using the         *
 *               non-universal rejection method for logconcave densities.    *
 *                                                                           *
 * REFERENCE:  - L. Devroye (1986): Non-Uniform Random Variate Generation,   *
 *               Springer Verlag, New York.                                  *
 *                                                                           *
 * Implemented by K. Lehner, 1990                                            *
 * Revised by F. Niederl, August 1992                                        *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define s    GEN.gen_param[0]
#define sm1  GEN.gen_param[1]

inline static void
powerexponential_epd_init( struct unur_gen *gen )
{
  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = 2;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  s = delta / 2.;
  sm1 = 1. - s;
  /* -X- end of setup code -X- */

} /* end of powerexponential_epd_init() */

double 
unur_stdgen_sample_powerexponential_epd( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double U,u1,V,X,y;
  double tau = 2./delta;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  while (1) {
    U = 2. * uniform() - 1.;                                  /* U(-1.0/1.0) */
    u1 = fabs(U);                                             /* u1=|u|      */
    V = uniform();                                            /* U(0/1)      */

    if (u1 <= sm1)
      /* Uniform hat-function for x <= (1-1/tau)   */
      X = u1;
    else {                       
      /* Exponential hat-function for x > (1-1/tau) */
      y = tau * (1. - u1);                                         /* U(0/1) */
      x = sm1 - s * log(y);
      V *= y;
    }
  } while (log(V) > -exp(log(X)*tau));               /* Acceptance/Rejection */
  
  /* Random sign */
  if (U < 0.)
    X = -X;

  /* -X- end of generator code -X- */

  return ((DISTR.n_params==1) ? X : theta + phi * X );

} /* end of unur_stdgen_sample_powerexponential_epd() */

/*---------------------------------------------------------------------------*/
#undef s
#undef sm1
/*---------------------------------------------------------------------------*/
