/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_weibull_gen.c                                              *
 *                                                                           *
 *   Special generators for Weibull distribution                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"

/*---------------------------------------------------------------------------*/
/* Prototypes for special generators                                         */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define c      (DISTR.params[0])   /* shape */
#define alpha  (DISTR.params[1])   /* scale */
#define zeta   (DISTR.params[2])   /* location */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_weibull_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Weibull distribution                */
     /* if gen == NULL then only check existance of variant.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case UNUR_STDGEN_INVERSION:   /* inversion method */
    if (par) PAR->is_inversion = TRUE;
    _unur_cstd_set_sampling_routine(par,gen,_unur_stdgen_sample_weibull_inv); 
    return UNUR_SUCCESS;

  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_weibull_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double _unur_stdgen_sample_weibull_inv( struct unur_gen *gen )
     /* Inversion method */
{
  /* -X- generator code -X- */
  double U,X;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INFINITY);

  /* sample from uniform random number generator */
  U = GEN->umin + uniform() * (GEN->umax-GEN->umin);

  /* transform to random variate */
  X = pow( -log(1.-U), 1./c );

  /* -X- end of generator code -X- */

  return ((DISTR.n_params==1) ? X : zeta + alpha * X );

} /* end of _unur_stdgen_sample_weibull_inv() */

/*---------------------------------------------------------------------------*/
