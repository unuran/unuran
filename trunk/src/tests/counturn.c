/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      counturn.c                                                   *
 *                                                                           *
 *   Count used uniform random numbers                                       *
 *                                                                           *
 *****************************************************************************
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
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
static char test_name[] = "Counting";

/*---------------------------------------------------------------------------*/
/* common variables                                                          */

static long urng_counter = 0;                /* count uniform random numbers */

/*---------------------------------------------------------------------------*/
/* wrapper for uniform random number generator that performs counting        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_FVOID 
/*---------------------------------------------------------------------------*/
static double (*urng_to_use)(void);          /* pointer to real uniform RNG  */
/*---------------------------------------------------------------------------*/

static double
_urng_with_counter(void)
     /*----------------------------------------------------------------------*/
     /* wrapper for uniform random number generator that performs counting   */
     /*----------------------------------------------------------------------*/
{
  ++urng_counter;
  return urng_to_use();
} /* end of urng_with_counter() */

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/
static double (*urng_to_use)(struct prng *); /* pointer to real uniform RNG  */
/*---------------------------------------------------------------------------*/

static double
_urng_with_counter( struct prng *gen )
     /*----------------------------------------------------------------------*/
     /* wrapper for uniform random number generator that performs counting   */
     /*----------------------------------------------------------------------*/
{
  ++urng_counter;
  return urng_to_use(gen);
} /* end of urng_with_counter() */

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/
static double (*urng_to_use)(void*);         /* pointer to real uniform RNG  */
/*---------------------------------------------------------------------------*/

static double
_urng_with_counter(void *params)
     /*----------------------------------------------------------------------*/
     /* wrapper for uniform random number generator that performs counting   */
     /*----------------------------------------------------------------------*/
{
  ++urng_counter;
  return urng_to_use(params);
} /* end of urng_with_counter() */

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_test_count_urn( struct unur_gen *gen, int samplesize, int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* count used uniform random numbers                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... sample size                                         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   total number of used uniform random numbers                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  long j;
  UNUR_URNG *urng_aux;

  /* check arguments */
  _unur_check_NULL(test_name,gen,-1);

  /* reset counter */
  urng_counter = 0;

  /* save auxilliary generator */
  urng_aux = gen->urng_aux;

  /* exchange pointer to uniform rng with counting wrapper */
#if UNUR_URNG_TYPE == UNUR_URNG_FVOID
  urng_to_use = gen->urng;
  unur_chg_urng(gen,_urng_with_counter);
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
  urng_to_use = gen->urng->get_next;
  gen->urng->get_next = _urng_with_counter;
  if (gen->urng_aux) gen->urng_aux = gen->urng;
#elif UNUR_URNG_TYPE == UNUR_URNG_GENERIC
  urng_to_use = gen->urng->getrand;
  gen->urng->getrand = _urng_with_counter;
  if (gen->urng_aux) gen->urng_aux = gen->urng;
#else
  /* no counter available */
  if (verbosity)
    fprintf(out,"\nCOUNT: ---  (cannot count urng for URNG type)\n");
  return -1;
#endif  /* UNUR_URNG_TYPE */

  /* run generator */
  switch (gen->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    for( j=0; j<samplesize; j++ )
      _unur_sample_discr(gen);
    break;
  case UNUR_METH_CONT:
    for( j=0; j<samplesize; j++ )
      _unur_sample_cont(gen);
    break;
  case UNUR_METH_VEC: 
    { /* we need an array for the vector */
      double *vec;
      int dim;
      dim = unur_get_dimension(gen);
      vec = _unur_malloc( dim * sizeof(double) );
      for( j=0; j<samplesize; j++ )
	_unur_sample_vec(gen,vec);
      free(vec);
    }
    break;
  default: /* unknown ! */
    _unur_error("Tests",UNUR_ERR_GENERIC,"method unknown!");
    return 0;
  }

  /* reset pointer to uniform rng */
#if UNUR_URNG_TYPE == UNUR_URNG_FVOID
  unur_chg_urng(gen,urng_to_use);
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
  gen->urng->get_next = urng_to_use;
#elif UNUR_URNG_TYPE == UNUR_URNG_GENERIC
  gen->urng->getrand = urng_to_use;
#endif  /* UNUR_URNG_TYPE */

  /* restore auxilliary generator */
  gen->urng_aux = urng_aux;

  /* print result */
  if (verbosity) {
    fprintf(out,"\nCOUNT: %g urng per generated number (total = %ld)\n",
	   ((double)urng_counter)/((double) samplesize),urng_counter);
  }

  /* return total number of used uniform random numbers */
  return urng_counter;

} /* end of unur_test_count_urn() */

/*---------------------------------------------------------------------------*/
