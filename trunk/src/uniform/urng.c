/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng.c                                                            *
 *                                                                           *
 *   unified interface for uniform random number generators                  *
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
#if UNUR_URNG_TYPE != UNUR_URNG_GENERIC
#else
/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include "unur_uniform.h"

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_new (double (*getrand)(void *params), double (*reset)(void *params), void *params)
     /*----------------------------------------------------------------------*/
     /* create a new object for uniform random number generator              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   debug ... debugging flag                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,par,UNUR_ERR_NULL );



double unur_urng_sample (UNUR_URNG *urng);
     /*----------------------------------------------------------------------*/
     /* set debugging flag for generator                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   debug ... debugging flag                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,par,UNUR_ERR_NULL );



int unur_urng_reset (UNUR_URNG *urng);
int unur_urng_free (UNUR_URNG *urng);



/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/

