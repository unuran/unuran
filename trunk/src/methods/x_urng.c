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

static UNUR_URNG *urng_default = NULL;
static UNUR_URNG *urng_aux_default = NULL;

/*---------------------------------------------------------------------------*/

UNUR_URNG *
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
#if UNUR_URNG_TYPE == UNUR_URNG_POINTER 
    urng_default = UNUR_URNG_DEFAULT;
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
    urng_default = prng_new(UNUR_URNG_DEFAULT);
    if( urng_default == NULL ) {
      /* some parameters invalid! */
      _unur_error("prng",UNUR_ERR_NULL,"Cannot set default URNG. EXIT !!!");
      /* we cannot recover from this error */
      exit(EXIT_FAILURE);
    }
#else
#error UNUR_URNG_TYPE not valid !!
#endif
  }

  /* return default generator */
  return (urng_default);
} /* end of unur_get_default_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_set_default_urng( UNUR_URNG *urng_new )
     /*----------------------------------------------------------------------*/
     /* set default uniform random number generator and return old one       */
     /*                                                                      */
     /* parameters: pointer to new default uniform random number generator   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to old  uniform random number generator                    */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng_old = urng_default;

  /* NULL pointer not allowed */
  _unur_check_NULL("URNG", urng_new, urng_default);

  urng_default = urng_new;     /* reset urng */

  /* return old default generator */
  return (urng_old);
} /* end of unur_set_default_urng() */

/*---------------------------------------------------------------------------*/

int
unur_use_urng_aux_default( UNUR_PAR *par )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator to default            */
     /* (initialize generator if necessary)                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  if (par->urng_aux == NULL)
    /* no auxilliary generator is required */
    return 0;

  /* default generator already running ? */
  if( urng_aux_default == NULL ) {
    /* have to initialize default generator first */
#if UNUR_URNG_TYPE == UNUR_URNG_POINTER 
    urng_aux_default = UNUR_URNG_AUX_DEFAULT;
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
    urng_aux_default = prng_new(UNUR_URNG_AUX_DEFAULT);
    if( urng_aux_default == NULL ) {
      /* some parameters invalid! */
      _unur_error("prng",UNUR_ERR_NULL,"Cannot set default aux URNG. EXIT !!!");
      /* we cannot recover from this error */
      return 0;
    }
#else
#error UNUR_URNG_TYPE not valid !!
#endif
  }

  /* set aux URNG */
  par->urng_aux = urng_aux_default;

  return 1;

} /* end of unur_use_urng_aux_default() */

/*---------------------------------------------------------------------------*/

int
unur_chgto_urng_aux_default( UNUR_GEN *gen )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator to default            */
     /* (initialize generator if necessary)                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  if (gen->urng_aux == NULL)
    /* no auxilliary generator is required */
    return 0;

  /* default generator already running ? */
  if( urng_aux_default == NULL ) {
    /* have to initialize default generator first */
#if UNUR_URNG_TYPE == UNUR_URNG_POINTER 
    urng_aux_default = UNUR_URNG_AUX_DEFAULT;
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
    urng_aux_default = prng_new(UNUR_URNG_AUX_DEFAULT);
    if( urng_aux_default == NULL ) {
      /* some parameters invalid! */
      _unur_error("prng",UNUR_ERR_NULL,"Cannot set default aux URNG. EXIT !!!");
      /* we cannot recover from this error */
      return 0;
    }
#else
#error UNUR_URNG_TYPE not valid !!
#endif
  }

  /* set aux URNG */
  gen->urng_aux = urng_aux_default;

  return 1;

} /* end of unur_chgto_urng_aux_default() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Set, get or change uniform RNG for generator                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_urng( struct unur_par *par, UNUR_URNG *urng )
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
  _unur_check_NULL( NULL,par,0 );
  _unur_check_NULL("URNG",urng,0);

  /* main generator */
  par->urng = urng;

  /* overwrite auxilliary generator */
  if (par->urng_aux) par->urng_aux = urng;

  return 1;
} /* end of unur_set_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
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

UNUR_URNG *
unur_chg_urng( struct unur_gen *gen, UNUR_URNG *urng )
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
  UNUR_URNG *urng_old;

  /* check arguments */
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng,NULL);

  urng_old = gen->urng;

  /* set pointer to main URNG */
  gen->urng = urng;

  /* also set pointer in auxilliary generator objects */
  if (gen->gen_aux)   unur_chg_urng(gen->gen_aux,urng);

  /* overwrite auxilliary URNG */
  if (gen->urng_aux) gen->urng_aux = urng;

  return urng_old;
} /* end of unur_chg_urng() */

/*****************************************************************************/

int
unur_set_urng_aux( struct unur_par *par, UNUR_URNG *urng_aux )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par_aux ... pointer to parameter for building generator object     */
     /*   urng    ... pointer to auxilliary uniform random number generator  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,par,0 );
  _unur_check_NULL("URNGaux",urng_aux,0);

  if (par->urng_aux == NULL)
    /* no auxilliary generator is required */
    return 0;

  par->urng_aux = urng_aux;

  return 1;
} /* end of unur_set_urng_aux() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_get_urng_aux( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get auxilliary uniform random number generator                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old auxilliary uniform RNG                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);

  return gen->urng_aux;
} /* end of unur_get_urng_aux() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_chg_urng_aux( struct unur_gen *gen, UNUR_URNG *urng_aux )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   urng_aux ... pointer to auxilliary uniform random number generator */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old auxilliary uniform RNG                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng_aux_old;

  /* check arguments */
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng_aux,NULL);

  if (gen->urng_aux == NULL) 
    /* no auxilliary generator is required */
    return NULL;

  urng_aux_old = gen->urng_aux;

  /* set pointer to main URNG */
  gen->urng_aux = urng_aux;

  /* also set pointer in auxilliary generator objects */
  if (gen->gen_aux)   unur_chg_urng_aux(gen->gen_aux,urng_aux);

  return urng_aux_old;
} /* end of unur_chg_urng_aux() */

/*---------------------------------------------------------------------------*/
