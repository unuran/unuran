/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_hypergeometric_gen.c                                       *
 *                                                                           *
 *   Special generators for Hypergeometric distribution                      *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.dstd        /* data for parameter object         */
#define GEN       gen->data.dstd        /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params   100     /** TODO:  maximal number of parameters for generator */

/* parameters */
#define N  params[0]
#define M  params[1]
#define n  params[2]

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int hypergeometric_xxxx_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_hypergeometric_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Hypergeometric distribution         */
     /* if gen == NULL then only check existance of variant.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case 1:  /** TODO: XXXX  method */
    //    _unur_dstd_set_sampling_routine( par,gen,_unur_stdgen_sample_hypergeometric_xxxx );
    //    return hypergeometric_xxxx_init( gen );
    return 0;

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_hypergeometric_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/** TODO: das ist der header aus WINRAND **/

/*****************************************************************************
 *                                                                           *
 *  Hypergeometric Distribution: XXXX method                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Hypergeometric             *
 *               distribution with parameters r (no. of failures given)      *
 *               and p (probability of success) valid for r > 0, 0 < p < 1.  *
 *               If G from Gamma(r) then K from Poiss(pG/(1-p)) is           *
 *               NB(r,p)-distributed.                                        *
 * REFERENCE:  - J.H. Ahrens, U. Dieter (1974): Computer methods for         *
 *               sampling from gamma, bet ??????                             *
 *               distributions, Computing 12, 223--246.                      *
 *                                                                           *
 * Implemented by  E. Stadlober,  September 1991                             *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

inline static int
hypergeometric_xxxx_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DSTD_GEN,0);

  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = MAX_gen_params;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */

  /** TODO: insert setup for generator **/

  /* -X- end of setup code -X- */

  return 1;

} /* end of hypergeometric_xxxx_init() */

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_hypergeometric_xxx( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double y;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DSTD_GEN,0);

  /** TODO: sample code **/

  /* -X- end of generator code -X- */

  return 0.;

} /* end of _unur_stdgen_sample_hypergeometric_xxxx() */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


