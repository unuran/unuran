/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_gsl.c                                                        *
 *                                                                           *
 *   routines to get new URNG object with sampling routine of type GSL.      *
 *   GSL (GNU Scientific Library), see http://www.gnu.org/software/gsl/.     *
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

#include <unur_source.h>
#include "urng.h"

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/
#ifdef UNURAN_HAS_GSL
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_gslptr_new( gsl_rng *gsl )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type GSL.                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gsl ... pointer to generator structure                             */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = unur_urng_new( (double(*)(void*)) gsl_rng_uniform_pos, gsl );
  unur_urng_set_delete(urng, (void(*)(void*)) gsl_rng_free);
  unur_urng_set_seed(urng, (void(*)(void*,unsigned long)) gsl_rng_set);
  return urng;
} /* end of unur_urng_gslptr_new() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_gsl_new( const gsl_rng_type *urngtype )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type PRNG.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prngstr ... string that describes generator                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = unur_urng_gslptr_new(gsl_rng_alloc(urngtype));
  unur_urng_seed(urng,gsl_rng_default_seed);
  return urng;
} /* end of unur_urng_gsl_new() */

/*---------------------------------------------------------------------------*/
#endif   /* #ifdef UNURAN_HAS_GSL */
/*---------------------------------------------------------------------------*/
#endif   /* #if UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/
