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

#include <unuran_tests.h>
#include <source_unuran.h>

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
static char test_name[] = "Counting";
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* wrapper for uniform random number generator that performs counting        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#if UNUR_URNG_INVOKE == UNUR_URNG_POINTER 
/*---------------------------------------------------------------------------*/
static long urng_counter = 0;                /* count uniform random numbers */
static double (*urng_to_use)(void);          /* pointer to real uniform rng  */
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
#elif UNUR_URNG_INVOKE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/
static long urng_counter = 0;                /* count uniform random numbers */
static double (*urng_to_use_next)(struct prng *); /* pointer to real URNG    */
/*---------------------------------------------------------------------------*/

static double
_urng_with_counter_next( struct prng *gen )
     /*----------------------------------------------------------------------*/
     /* wrapper for uniform random number generator that performs counting   */
     /*----------------------------------------------------------------------*/
{
  ++urng_counter;
  return urng_to_use_next(gen);
} /* end of urng_with_counter() */

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_INVOKE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_INVOKE */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_test_count_urn( struct unur_gen *gen, int samplesize )
     /*----------------------------------------------------------------------*/
     /* count used uniform random numbers                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   total number of used uniform random numbers                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  long j;

  /* check arguments */
  _unur_check_NULL(test_name,gen,-1);

  /* reset counter */
  urng_counter = 0;

  /* exchange pointer to uniform rng with counting wrapper */
#if UNUR_URNG_INVOKE == UNUR_URNG_POINTER
  urng_to_use = gen->urng;
  unur_chg_urng(gen,_urng_with_counter);
#elif UNUR_URNG_INVOKE == UNUR_URNG_PRNG
  urng_to_use_next = gen->urng->get_next;
  gen->urng->get_next = _urng_with_counter_next;
#endif

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
#if UNUR_URNG_INVOKE == UNUR_URNG_POINTER
  unur_chg_urng(gen,urng_to_use);
#elif UNUR_URNG_INVOKE == UNUR_URNG_PRNG
  gen->urng->get_next = urng_to_use_next;
#endif

  /* print result */
  printf("\nCOUNT: %g urng per generated number (total = %ld)\n",
	 ((double)urng_counter)/((double) samplesize),urng_counter);

  /* return total number of used uniform random numbers */
  return urng_counter;
} /* end of unur_test_count_urn() */
