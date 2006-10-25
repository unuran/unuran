/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_burr_gen.c                                                 *
 *                                                                           *
 *   Special generators for Burr family of distributions                     *
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
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/
/* Prototypes for special generators                                         */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define burr_type (gen->distr->id)
#define k         (DISTR.params[1])
#define c         (DISTR.params[2])

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_burr_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Burr family of distributions        */
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
    if (par->distr->id == UNUR_DISTR_BURR_XI) {
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (gen) GEN->is_inversion = TRUE;
    _unur_cstd_set_sampling_routine(par,gen,_unur_stdgen_sample_burr_inv); 
    return UNUR_SUCCESS;

  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_burr_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double _unur_stdgen_sample_burr_inv( struct unur_gen *gen )
     /* Inversion method                                                     */
{
  /* -X- generator code -X- */
  double U, Y;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  /* sample from uniform random number generator */
  while ((U = GEN->umin + uniform() * (GEN->umax - GEN->umin)) == 0.);

  /* transform to random variate */
  switch (burr_type) {

  case UNUR_DISTR_BURR_I:
    return U;

  case UNUR_DISTR_BURR_II:
    Y = exp( -log(U)/k );  /* U^(-1/k) */
    return ( -log( Y - 1. ) );

  case UNUR_DISTR_BURR_III:
    Y = exp( -log(U)/k );  /* U^(-1/k) */
    return ( exp( -log( Y - 1. )/c ) );

  case UNUR_DISTR_BURR_IV:
    Y = exp( -log(U)/k );   /* U^(-1/k) */
    Y = exp( c * log( Y - 1. )) + 1.;
    return (c/Y);

  case UNUR_DISTR_BURR_V:
    Y = exp( -log(U)/k );   /* U^(-1/k) */
    return atan( -log( (Y - 1.) / c ) );

  case UNUR_DISTR_BURR_VI:
    Y = exp( -log(U)/k );   /* U^(-1/k) */
    Y = -log( (Y - 1.) / c)/k;
    return log( Y + sqrt(Y * Y +1.));

  case UNUR_DISTR_BURR_VII:
    Y = exp( log(U)/k );    /* U^(1/k) */
    return ( log(2. * Y / (2. - 2.*Y)) / 2. );

  case UNUR_DISTR_BURR_VIII:
    Y = exp( log(U)/k );    /* U^(1/k) */
    return ( log( tan( Y * M_PI/2. ) ) );

  case UNUR_DISTR_BURR_IX:
    Y = 1. + 2. * U / (c * (1.-U));
    return log( exp( log(Y) / k) - 1. );

  case UNUR_DISTR_BURR_X:
  Y = exp( log(U)/k );   /* U^(1/k) */
    return ( sqrt( -log( 1. - Y ) ) );

  case UNUR_DISTR_BURR_XII:
    Y = exp( -log(U)/k );   /* U^(-1/k) */
    return ( exp( log( Y - 1.) / c) );

  case UNUR_DISTR_BURR_XI:
  default:
    _unur_error(NULL,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }

  /* -X- end of generator code -X- */

} /* end of _unur_stdgen_sample_burr_inv() */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#undef burr_type
#undef k
#undef c
/*---------------------------------------------------------------------------*/
