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
