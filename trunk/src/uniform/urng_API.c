/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_API.c                                                        *
 *                                                                           *
 *   Unified interface for UNURAN uniform random number generators.          *
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
#include "unur_uniform.h"
#include "urng.h"


/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Handle and sample from URNG object                                     **/
/**                                                                         **/
/*****************************************************************************/

double
unur_urng_sample (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Sample from URNG object.                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   delete ... function for destroying generator object                */
     /*                                                                      */
     /* return:                                                              */
     /*   uniform random number                                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  return _unur_call_urng(urng);
} /* end of unur_urng_sample() */ 

/*---------------------------------------------------------------------------*/


int
unur_urng_reset (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Reset URNG object.                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC

  /* check whether we can reset the URNG object */
  if (urng->reset == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"reset");
    return UNUR_ERR_URNG_MISS;
  }

  /* reset object */
  urng->reset (urng->params);
  /** TODO: is there any chance to have an exitcode from this call? */

  return UNUR_SUCCESS;

#else

  return _unur_call_reset(urng);

#endif

} /* end of unur_urng_reset() */ 

/*---------------------------------------------------------------------------*/








/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Create new URNG object                                                 **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_new( double (*sampleunif)(void *p), void *params )
     /*----------------------------------------------------------------------*/
     /* create a new URNG object for uniform random number generator.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   sampleunif ... pointer to sampling routine                         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to URNG object                                             */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = NULL;         /* pointer to URNG object */

  /* check arguments */
  _unur_check_NULL( "URNG", sampleunif, NULL );

  /* allocate memory for URNG object */
  urng = _unur_xmalloc( sizeof(struct unur_urng_generic) );

  /* copy parameters into object */
  urng->sampleunif = sampleunif;
  urng->params  = params;

  /* initialize optional functions (set to not available) */
  urng->delete = NULL;
  urng->reset = NULL;
  urng->nextsub = NULL;
  urng->resetsub = NULL;
  urng->anti = NULL;

  /* return object */
  return urng;
} /* end of unur_urng_new() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_anti( UNUR_URNG *urng, int (*anti)(void *p, int anti) )
     /*----------------------------------------------------------------------*/
     /* Set function to switch to antithetic random numbers.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   anti   ... function for switching                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );

  urng->anti = anti;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_anti() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_reset( UNUR_URNG *urng, int (*reset)(void *p) )
     /*----------------------------------------------------------------------*/
     /* Set function for reseting URNG.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   reset  ... function for reseting generator object                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );

  urng->reset = reset;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_reset() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_nextsub( UNUR_URNG *urng, int (*nextsub)(void *p) )
     /*----------------------------------------------------------------------*/
     /* Set function for jumping to start of the next substream of URNG.     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng    ... pointer to URNG object                                 */
     /*   nextsub ... function for jumping to next substream                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );

  urng->nextsub = nextsub;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_nextsub() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_resetsub( UNUR_URNG *urng, int (*resetsub)(void *p) )
     /*----------------------------------------------------------------------*/
     /* Set function for jumping to start of the current substream of URNG.  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng     ... pointer to URNG object                                */
     /*   resetsub ... function for reseting current substream               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );

  urng->resetsub = resetsub;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_resetsub() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_delete( UNUR_URNG *urng, void (*delete)(void *p) )
     /*----------------------------------------------------------------------*/
     /* Set function for destroying URNG.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   delete ... function for destroying generator object                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );

  urng->delete  = delete;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_delete() */

/*---------------------------------------------------------------------------*/

int
unur_urng_anti (UNUR_URNG *urng, int anti)
     /*----------------------------------------------------------------------*/
     /* Set antithetic flag for URNG.                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   anti   ... antithetic flag                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  /* check whether we can set the antithetic flag */
  if (urng->anti == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"antithetic flag");
    return UNUR_ERR_URNG_MISS;
  }

  /* set flag */
  urng->anti (urng->params,anti);
  /** TODO: is there any chance to have an exitcode from this call? */

  return UNUR_SUCCESS;
} /* end of unur_urng_anti() */ 

/*---------------------------------------------------------------------------*/

int
unur_urng_nextsub (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Jump to start of next substream in URNG.                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  /* check whether we can reset the URNG object */
  if (urng->nextsub == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"next substream");
    return UNUR_ERR_URNG_MISS;
  }

  /* jump to next substream */
  urng->nextsub (urng->params);
  /** TODO: is there any chance to have an exitcode from this call? */

  return UNUR_SUCCESS;
} /* end of unur_urng_nextsub() */ 

/*---------------------------------------------------------------------------*/

int
unur_urng_resetsub (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Reset current substream of the URNG object.                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  /* check whether we can reset the URNG object */
  if (urng->resetsub == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"reset substream");
    return UNUR_ERR_URNG_MISS;
  }

  /* reset substream */
  urng->resetsub (urng->params);
  /** TODO: is there any chance to have an exitcode from this call? */

  return UNUR_SUCCESS;
} /* end of unur_urng_resetsub() */ 

/*---------------------------------------------------------------------------*/

int
unur_urng_free (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Destroy URNG object.                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );

  if (urng->delete != NULL) urng->delete (urng->params);
  free (urng);
  urng = NULL;

  return UNUR_SUCCESS;
} /* end of unur_urng_free() */ 

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Handle URNG object in a generator object                               **/
/**                                                                         **/
/*****************************************************************************/

int
unur_gen_anti (UNUR_GEN *gen, int anti)
     /*----------------------------------------------------------------------*/
     /* Set antithetic flag of uniform generator in generator object.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*   anti ... antithetic flag                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );

  return unur_urng_anti(gen->urng, anti);
} /* end of unur_gen_reset() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_reset (UNUR_GEN *gen)
     /*----------------------------------------------------------------------*/
     /* Reset uniform generator in generator object.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_reset(gen->urng);
} /* end of unur_gen_reset() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_nextsub (UNUR_GEN *gen)
     /*----------------------------------------------------------------------*/
     /* Jump uniform generator in generator object to start of next          */
     /* substream.                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );

  return unur_urng_nextsub(gen->urng);
} /* end of unur_gen_nextsub() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_resetsub (UNUR_GEN *gen)
     /*----------------------------------------------------------------------*/
     /* Reset current substream of uniform generator in generator object.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );

  return unur_urng_resetsub(gen->urng);
} /* end of unur_gen_resetsub() */ 

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/

