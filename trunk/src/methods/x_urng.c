/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_urng.c                                                          *
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

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Default uniform random number generator                                **/
/**                                                                         **/
/*****************************************************************************/

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
#if UNUR_URNG_INVOKE == UNUR_URNG_POINTER 
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

/*****************************************************************************/
/**                                                                         **/
/**  Set, get or change uniform RNG for generator                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_urng( struct unur_par *par, UNUR_URNG_TYPE urng )
     /*----------------------------------------------------------------------*/
     /* set uniform random number generator                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   urng    ... pointer to uniform random number generator             */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);
  CHECK_NULL(urng,0);

  par->urng = urng;

  return 1;
} /* end of unur_set_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG_TYPE
unur_get_urng( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get uniform random number generator                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old uniform RNG                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);

  return gen->urng;

} /* end of unur_get_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG_TYPE
unur_chg_urng( struct unur_gen *gen, UNUR_URNG_TYPE urng )
     /*----------------------------------------------------------------------*/
     /* set uniform random number generator                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   urng    ... pointer to uniform random number generator             */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old uniform RNG                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG_TYPE urng_old;

  /* check arguments */
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng,NULL);

  urng_old = gen->urng;

  gen->urng = urng;

  return urng_old;
} /* end of unur_chg_urng() */

/*---------------------------------------------------------------------------*/


