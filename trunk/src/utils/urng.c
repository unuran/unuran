/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      urng.c                                                       *
 *                                                                           *
 *   routines for invoking uniform random number generator                   *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Thu Jul 22 16:07:30 CEST 1999                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_urng.h>
#include <unur_errno.h>

/*---------------------------------------------------------------------------*/

/* pointer to default uniform random number generator */
static UNUR_URNG_TYPE urng_default = NULL;

/*---------------------------------------------------------------------------*/

UNUR_URNG_TYPE
unur_get_default_urng( void )
     /*----------------------------------------------------------------------*/
     /* return default uniform random number generator                       */
     /* (initialize generator if necessary)                                  */
     /*                                                                      */
     /* parameters: none                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to default generator                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* default generator already running ? */
  if( urng_default == NULL ) {
    /* have to initialize default generator first */
#if UNUR_URNG_INVOKE == UNUR_URNG_LINKED
    ;    /* this function is not required for this compiler switch. */
#elif UNUR_URNG_INVOKE == UNUR_URNG_POINTER 
    urng_default = UNUR_URNG_DEFAULT;
#elif UNUR_URNG_INVOKE == UNUR_URNG_PRNG
    urng_default = prng_new(UNUR_URNG_DEFAULT);
    if( urng_default == NULL )
      /* some parameters invalid! */
      _unur_error("prng",UNUR_ERR_NULL,"Cannot set default URNG");
#else
#error UNUR_URNG_INVOKE not valid !!
#endif
  }

  /* return default generator */
  return (urng_default);
} /* end of unur_get_default_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG_TYPE
unur_set_default_urng( UNUR_URNG_TYPE urng_new )
     /*----------------------------------------------------------------------*/
     /* set default uniform random number generator and return old one       */
     /*                                                                      */
     /* parameters: pointer to new default uniform random number generator   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to old  uniform random number generator                    */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG_TYPE urng_old = urng_default;

  /* NULL pointer not allowed */
  if (urng_new == NULL) {
    _unur_warning(NULL,UNUR_ERR_NULL,"invalid URNG.");
    return urng_default;
  }

  urng_default = urng_new;     /* reset urng */

  /* return old default generator */
  return (urng_old);
} /* end of unur_set_default_urng() */

/*---------------------------------------------------------------------------*/
