/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_geometric_gen.c                                            *
 *                                                                           *
 *   Special generators for Geometric distribution                           *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.dstd        /* data for parameter object         */
#define GEN       gen->data.dstd        /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

/* parameters */
#define p  (DISTR.params[0])    /* shape */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_geometric_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Geometric distribution              */
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
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case UNUR_STDGEN_INVERSION:   /* inversion method */
    if (par) PAR.is_inversion = TRUE;
    _unur_dstd_set_sampling_routine( par,gen,_unur_stdgen_sample_geometric_inv );
    return 1;

  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_geometric_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Geometric Distribution: Inversion                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Geometric distribution     *
 *               with parameter 0 < p < 1.                                   *
 *                                                                           *
 * Implemented by R.Kremer 1990, revised by P.Busswald, July 1992            *
 *****************************************************************************
 *                                                                           *
 *  On generating random numbers of a discrete distribution by Inversion     *
 *  normally sequential search is necessary, but in the case of the          *
 *  Geometric distribution a direct transformation is possible because of    *
 *  the special parallel to the continuous Exponential distribution Exp(t):  *
 *    X - Exp(t): G(x)=1-exp(-tx)                                            *
 *        Geo(p): pk=G(k+1)-G(k)=exp(-tk)*(1-exp(-t))                        *
 *                p=1-exp(-t)                                                *
 *  A random number of the Geometric distribution Geo(p) is obtained by      *
 *  k = (int)x, where x is from Exp(t) with parameter t = -log(1-p).         *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_geometric_inv( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double U;
  int K;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DSTD_GEN,0);

  /* sample from uniform random number generator */
/*    while ((U = GEN.umin + uniform() * (GEN.umax-GEN.umin)) == 0.); */
  while ((U = uniform()) == 0.);

  /* transform to random variate */
  K = (int) (log(U) / log(1.-p));

  /* -X- end of generator code -X- */

  return K;
  
} /* end of _unur_stdgen_sample_geometric_inv() */

/*---------------------------------------------------------------------------*/
