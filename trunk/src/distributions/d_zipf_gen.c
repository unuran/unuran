/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_zipf_gen.c                                                 *
 *                                                                           *
 *   Special generators for Zipf (or Zeta) distribution                      *
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
#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static void zipf_zet_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.dstd        /* data for parameter object         */
#define GEN       gen->data.dstd        /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

/* parameters */
#define rho  (DISTR.params[0])    /* shape */
#define tau  (DISTR.params[1])

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_zipf_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Zipf distribution                    */
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
  CHECK_NULL(par,0.);
  COOKIE_CHECK(par,CK_DSTD_PAR,0.);

  switch (par->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* Acceptance Rejection */
    _unur_dstd_set_sampling_routine( par,gen,unur_stdgen_sample_zipf_zet );
    zipf_zet_init( gen );
    return 1;

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(par->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_zipf_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Zipf Distribution: Acceptance Rejection                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Zipf distribution with     *    
 *               parameters rho > 0 and tau >= 0.                            *
 *                                                                           *
 * REFERENCE : - J. Dagpunar (1988): Principles of Random Variate Generation,*
 *               Clarendon Press, Oxford.                                    *
 *                                                                           *
 * Implemented by P. Busswald, September 1992                                *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#if 0
#define b       (GEN.gen_param[0])
#define vm      (GEN.gen_param[1])
#define vp      (GEN.gen_param[2])
#define vd      (GEN.gen_param[3])
#endif
/*---------------------------------------------------------------------------*/

inline static void
zipf_zet_init( struct unur_gen *gen )
{
#if 0
  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = 4;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  if (nu < 1.) {
    _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
    return;
  }

  if (nu == 1.)
    /* no setup step required */
    return;

  /* else nu > 1 */
  b = sqrt(nu - 1.);
  vm = - 0.6065306597 * (1. - 0.25 / (b * b + 1.));
  vm = (-b > vm) ? -b : vm;
  vp = 0.6065306597 * (0.7071067812 + b) / (0.5 + b);
  vd = vp - vm;
  /* -X- end of setup code -X- */
#endif
  ;
} /* end of zipf_chru_init() */

#define ro  (DISTR.params[0])    /* shape */
#define pk  (DISTR.params[1])

int
unur_stdgen_sample_zipf_zet( struct unur_gen *gen )
{
  /* -X- generator code -X- */

  static double c,d,ro_prev = -1.0,pk_prev = -1.0, maxlongint=2147483647.5;
         /*max. long-Zahl ist eigentl. maxlongint-1.5*/
  double u,v,e,x;
  long k;

  if (ro != ro_prev || pk != pk_prev)                    /* Set-up */
         {
          ro_prev = ro;
          pk_prev = pk;
          if (ro<pk)
        {
         c = pk-0.5;
         d = 0;
        }
          else
        {
         c = ro-0.5;
         d = (1.0+ro)*log((1.0+pk)/(1.0+ro));
        }
         }
  do
          {
                 do
         {
                u=uniform();
                v=uniform();
                x = (c+0.5)*exp(-log(u)/ro) - c;
                }
                 while (x<=0.5 || x>=maxlongint);
                 k = (long int) (x+0.5);
                 e = -log(v);
                }
  while ( e < (1.0+ro)*log((k+pk)/(x+c)) - d );
  return(k);


  /* -X- end of generator code -X- */
  

} /* end of unur_stdgen_sample_zipf_zet() */

/*---------------------------------------------------------------------------*/
#undef b 
#undef vm
#undef vp
#undef vd
/*---------------------------------------------------------------------------*/
