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

#define MAX_gen_params   1     /* maximal number of parameters for generator */

/* parameters */
#define p  (DISTR.params[0])
#define r  (DISTR.params[1])

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int negativebinomial_nbp_init( struct unur_gen *gen );

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
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* Compound method */
    _unur_dstd_set_sampling_routine( par,gen,_unur_stdgen_sample_negativebinomial_nbp );
    return negativebinomial_nbp_init( gen );

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
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
#define x   (GEN.gen_param[0])

#define GAMMA    gen->gen_aux    /* pointer to gamma variate generator       */
#define POISSON  gen->gen_aux_2  /* pointer to poisson variate generator     */
/*---------------------------------------------------------------------------*/

inline static int
negativebinomial_nbp_init( struct unur_gen *gen )
{
  double gamma_param;   
  double poisson_param;   

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DSTD_GEN,0);

  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = MAX_gen_params;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  x = p /(1. - p);

  /* make a gamma variate generator (use default special generator) */
  gamma_param = r;   /* shape parameter for gamma distribution */
  if (GAMMA==NULL) {
    struct unur_par *par = unur_cstd_new( unur_distr_gamma(&gamma_param,1) );
    GAMMA = (par) ? _unur_init(par) : NULL;
    _unur_check_NULL( NULL,GAMMA,0 );
    /* need same uniform random number generator as negative binomial generator */
    GAMMA->urng = gen->urng;
    /* copy debugging flags */
    GAMMA->debug = gen->debug;
  }
  else  /* re-init mode --> change shape parameter */
    unur_cstd_chg_pdfparams(GAMMA,&gamma_param,1);

  /* make a poisson variate generator (use default special generator) */
  if (POISSON==NULL) {
    poisson_param = 1.;   /* shape parameter for poisson distribution (use a dummy value yet) */
    POISSON = _unur_dstd_init( unur_dstd_new( unur_distr_poisson(&poisson_param,1) ) );
    _unur_check_NULL( NULL,POISSON,0 );
    /* need same uniform random number generator as negative binomial generator */
    POISSON->urng = gen->urng;
    /* copy debugging flags */
    POISSON->debug = gen->debug;
  }
  /* else we are in the re-init mode 
     --> there is no necessity to make the generator object again */

  /* -X- end of setup code -X- */

  return 1;

} /* end of poisson_pdac_init() */

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_negativebinomial_nbp( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double y;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DSTD_GEN,0);

  /* sample from gamma distribution */
  y = x * _unur_sample_cont(GAMMA);

  /* sample from poisson distribution */
  unur_dstd_chg_pdfparams(POISSON,&y,1);
  return _unur_sample_discr(POISSON);

  /* -X- end of generator code -X- */

} /* end of _unur_stdgen_sample_negativebinomial_nbp() */

/*---------------------------------------------------------------------------*/
#undef x

#undef GAMMA
#undef POISSON
/*---------------------------------------------------------------------------*/


