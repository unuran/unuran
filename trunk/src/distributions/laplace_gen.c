/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      laplace_gen.c                                                *
 *                                                                           *
 *   Special generators for Laplace distribution                             *
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
/* Prototypes for special generators                                         */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.cstd        /* data for parameter object         */
#define GEN       gen->data.cstd        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

#define theta  (DISTR.params[0])   /* location */
#define phi    (DISTR.params[1])   /* scale */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_laplace_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Laplace distribution                */
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

  case 0: /* Default */
  case UNUR_STDGEN_INVERSION:   /* inversion method */
    PAR.is_inversion = TRUE;
    _unur_cstd_set_sampling_routine(par,gen,unur_stdgen_sample_laplace_inv); 
    return 1;

  default:
    /* no such generator */
    _unur_warning(par->genid,UNUR_ERR_DISTR_GEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_laplace_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_laplace_inv( struct unur_gen *gen )
     /* Inversion method */
{
  /* -X- generator code -X- */
  double U, X;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  /* sample a uniform random number != 0 */
  while ((U = 2. * uniform()) == 0. ) ;

  X = (U>1.) ? -log(2.-U) : log(U);

  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : theta + phi * X );

} /* end of unur_stdgen_sample_laplace_inv() */

/*---------------------------------------------------------------------------*/
