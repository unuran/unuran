/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      get.h                                                        *
 *                                                                           *
 *   get parameters from (existing) generators                               *
 *                                                                           *
 *   PARAMETER: struct unur_gen *                                            *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Fri Dec 17 09:55:33 CET 1999                         *
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

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for the distribution and its p.d.f.                         **/
/**                                                                         **/
/*****************************************************************************/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for generators of discrete distributions                    **/
/**                                                                         **/
/*****************************************************************************/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for generators of continuous distributions                  **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
unur_get_n_intervals( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get number of intervals/segments                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   number of intervals/segments                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,-1);

  switch (gen->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(gen,CK_AROU_GEN,-1);
    return gen->data.arou.n_segs;

  case UNUR_METH_TABL:
    COOKIE_CHECK(gen,CK_TABL_GEN,-1);
    return gen->data.tabl.n_ivs;

  case UNUR_METH_TDR:
    COOKIE_CHECK(gen,CK_TDR_GEN,-1);
    return gen->data.tdr.n_ivs;

  default:
    /** TODO: error message **/
    return -1;
  }
  
} /* end of unur_get_n_intervals() */

/*---------------------------------------------------------------------------*/

int 
unur_get_max_intervals( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get maximal number of intervals/segments                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   maximal number of intervals/segments                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,-1);

  switch (gen->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(gen,CK_AROU_GEN,-1);
    return gen->data.arou.max_segs;

  case UNUR_METH_TABL:
    COOKIE_CHECK(gen,CK_TABL_GEN,-1);
    return gen->data.tabl.max_ivs;

  case UNUR_METH_TDR:
    COOKIE_CHECK(gen,CK_TDR_GEN,-1);
    return gen->data.tdr.max_ivs;

  default:
    /** TODO: error message **/
    return -1;
  }
  
} /* end of unur_get_max_intervals() */

/*---------------------------------------------------------------------------*/

double
unur_get_shratio( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get ratio A(squeeze) / A(hat)                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   ratio                                                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1.                                                         */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,-1.);

  switch (gen->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(gen,CK_AROU_GEN,-1.);
    return gen->data.arou.Asqueeze / gen->data.arou.Atotal;   /** TODO check for overflow ?? **/

  case UNUR_METH_TABL:
    COOKIE_CHECK(gen,CK_TABL_GEN,-1.);
    return gen->data.tabl.Asqueeze / gen->data.tabl.Atotal;   /** TODO check for overflow ?? **/

  case UNUR_METH_TDR:
    COOKIE_CHECK(gen,CK_TDR_GEN,-1.);
    return gen->data.tdr.Asqueeze / gen->data.tdr.Atotal;   /** TODO check for overflow ?? **/

  default:
    /** TODO: error message **/
    return -1.;
  }
  
} /* end of unur_get_shratio() */

/*---------------------------------------------------------------------------*/

int
unur_get_dimension( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* get dimension of generator for multivariate distribution                  */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*                                                                           */
/* return:                                                                   */
/*   dimension of distribution                                               */
/*---------------------------------------------------------------------------*/
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
    return 1;
  }

} /* end of unur_get_dimension() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for all generators                                          **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

UNUR_URNG_TYPE
unur_get_urng( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* set uniform random number generator                                       */
/*                                                                           */
/* parameters:                                                               */
/*   gen     ... pointer to generator object                                 */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);

  return gen->urng;

} /* end of unur_get_urng() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Auxilliary routines                                                    **/
/**                                                                         **/
/*****************************************************************************/

