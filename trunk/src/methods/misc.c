/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      misc.c                                                       *
 *                                                                           *
 *   miscelleanous routines                                                  *
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

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Set debuging flags                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_debug( struct unur_par *par, unsigned long debug )
     /*----------------------------------------------------------------------*/
     /* set debugging flag for generator                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   debug ... debugging flag                                           */
     /*----------------------------------------------------------------------*/
{
  CHECK_NULL(par,0);

#if UNUR_DEBUG & UNUR_DB_INFO
  par->debug = debug;
  return 1;
#else
  _unur_warning("DEBUG",UNUR_ERR_SET,"Debugging not available. Recompile library.");
  return 0;
#endif

} /* end of unur_set_debug() */
  
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

/*****************************************************************************/
/**                                                                         **/
/**  Get data about generator object                                        **/
/**                                                                         **/
/*****************************************************************************/

int
unur_get_dimension( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get dimension of generator for multivariate distribution             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   dimension of distribution                                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);

  switch (gen->method & UNUR_MASK_METHOD) {
  case UNUR_METH_RECT:
    COOKIE_CHECK(gen,CK_RECT_GEN,0);
    return gen->data.rect.dim;
    break;
  default:
    /* method unknown assume dim = 1 */
    /** TODO: distinguish between univariate --> dim = 1 and
	unknown multivariate **/
    return 1.;
  }

} /* end of unur_get_dimension() */

/*---------------------------------------------------------------------------*/



