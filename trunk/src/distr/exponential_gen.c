/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      exponential.c                                                *
 *                                                                           *
 *   Special generators for Exponential distribution                         *
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
#include <unur_distr.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
/* Prototypes for special generators                                         */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define GEN     gen->data.cstd
#define uniform()  (_unur_call_urng(gen))

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  get special sampling routine for distribution                          **/
/**                                                                         **/
/*****************************************************************************/

_UNUR_SAMPLING_ROUTINE_CONT *
_unur_stdgen_exponential_get_routine(unsigned variant)
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
  case UNUR_STDGEN_INVERSION:
    return unur_stdgen_sample_exponential_inv; /* inversion method */
  default:
    return NULL;
  }

} /* end of _unur_stdgen_exponential_get_routine() */

/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO

const char *
_unur_stdgen_exponential_routinename(void *routine)
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

  routinename( unur_stdgen_sample_exponential_inv );

  return NULL;

} /* end of _unur_stdgen_exponential_routinename() */

#endif

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Wrapper for special generators (WinRand)                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_exponential_inv( struct unur_gen *gen )
     /* Inversion method                                                     */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return ( -GEN.pdf_param[0] * log(1. - uniform()) );

} /* end of unur_stdgen_sample_exponential_inv() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

