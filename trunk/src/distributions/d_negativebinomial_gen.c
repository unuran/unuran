/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_negativebinomial_gen.c                                     *
 *                                                                           *
 *   Special generators for Negative Binomial distribution                   *
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

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

/* parameters */
#define p  (DISTR.params[0])
#define r  (DISTR.params[1])

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static void negativebinomial_nbp_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_negativebinomial_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Negative Binomial distribution      */
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
  /* check arguments */
  CHECK_NULL(par,0);  COOKIE_CHECK(par,CK_DSTD_PAR,0);

  switch (par->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* Compound method */
    _unur_dstd_set_sampling_routine( par,gen,unur_stdgen_sample_negativebinomial_nbp );
    negativebinomial_nbp_init( gen );
    return 1;

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(par->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_negativebinomial_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Negative Binomial Distribution: Compound method                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Negative Binomial          *
 *               distribution with parameters r (no. of failures given)      *
 *               and p (probability of success) valid for r > 0, 0 < p < 1.  *
 *               If G from Gamma(r) then K from Poiss(pG/(1-p)) is           *
 *               NB(r,p)-distributed.                                        *
 * REFERENCE:  - J.H. Ahrens, U. Dieter (1974): Computer methods for         *
 *               sampling from gamma, beta, Poisson and binomial             *
 *               distributions, Computing 12, 223--246.                      *
 *                                                                           *
 * Implemented by  E. Stadlober,  September 1991                             *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define GAMMA  GEN.gen_aux    /* pointer to gamma variate generator          */
/*---------------------------------------------------------------------------*/

inline static void
negativebinomial_nbp_init( struct unur_gen *gen );
{
  struct unur_par *par;

  /* check arguments */
  CHECK_NULL(gen,/*void*/); COOKIE_CHECK(gen,CK_DSTD_GEN,/*void*/);

  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = MAX_gen_params;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
    GEN.n_gen_iparam = MAX_gen_iparams;
    GEN.gen_iparam = _unur_malloc(GEN.n_gen_param * sizeof(int));
  }

  /* -X- setup code -X- */

  /* make a normal variate generator (use default special generator) */
  par = unur_cstd_new( unur_distr_gamma(NULL,0) );
  GAMMA = unur_init( par );
  /* need same uniform random number generator as slash generator */
  GAMMA->urng = gen->urng;

  /* -X- end of setup code -X- */

} /* end of poisson_pdac_init() */

/*---------------------------------------------------------------------------*/

int
unur_stdgen_sample_negativebinomial_nbp( struct unur_gen *gen )
{
  /* -X- generator code -X- */
 double y;
 static double x,p1 = -1.0;

  /* check arguments */
  CHECK_NULL(gen,0.); COOKIE_CHECK(gen,CK_DSTD_GEN,0.);

 if (p1 != p)
    {
     x = p /(1.0 - p);
     p1 = p;
    }
 y = x * Gamma(r,seed);
 return(Poisson(y,seed));

  /* -X- end of generator code -X- */

} /* end of unur_stdgen_sample_negativebinomial_nbp() */

/*---------------------------------------------------------------------------*/
