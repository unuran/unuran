/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_uniform_gen.c                                              *
 *                                                                           *
 *   Special generators for Uniform distribution                             *
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
/* Prototypes for special generators                                         */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.cstd        /* data for parameter object         */
#define GEN       gen->data.cstd        /* data for generator object         */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define a  (DISTR.params[0])    /* location */
#define b  (DISTR.params[1])    /* location */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_uniform_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Uniform distribution                */
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
    _unur_cstd_set_sampling_routine(par,gen,_unur_stdgen_sample_uniform_inv); 
    return 1;

  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_uniform_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double _unur_stdgen_sample_uniform_inv( struct unur_gen *gen )
     /* Inversion method */
{
  /* -X- generator code -X- */
  double X;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  /* sample from uniform random number generator */
  X = GEN.umin + uniform() * (GEN.umax-GEN.umin);

  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : a + X * (b - a));

} /* end of _unur_stdgen_sample_uniform_inv() */

/*---------------------------------------------------------------------------*/
