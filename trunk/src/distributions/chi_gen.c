/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      chis_gen.c                                                   *
 *                                                                           *
 *   Special generators for Chi distribution                                 *
 *   (Do not confuse with Chi^2 distribution!)                               *
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

inline static void chi_chru_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.cstd        /* data for parameter object         */
#define GEN       gen->data.cstd        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

#define nu (DISTR.params[0])    /* shape (degrees of freedom) */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_chi_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for chi distribution                    */
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
  COOKIE_CHECK(par,CK_CSTD_PAR,0.);

  switch (par->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* Ratio of Uniforms with shift */
    if (par->distr->data.cont.params[0] < 1.) {    /* nu < 1 */
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return 0;
    }
    /* nu >= 1 !!!! */
    _unur_cstd_set_sampling_routine( par,gen,unur_stdgen_sample_chi_chru );
    chi_chru_init( gen );
    return 1;

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(par->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_chi_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Chi Distribution: Ratio of Uniforms with shift                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Chi distribution with      *    
 *               parameters nu >= 1.                                         *
 *                                                                           *
 * REFERENCE : - J.F. Monahan (1987): An algorithm for generating chi random *
 *               variables, ACM Trans. Math. Software 13, 168-172.           *
 *                                                                           *
 * Implemented by R. Kremer, 1990                                            *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define b       (GEN.gen_param[0])
#define vm      (GEN.gen_param[1])
#define vp      (GEN.gen_param[2])
#define vd      (GEN.gen_param[3])
/*---------------------------------------------------------------------------*/

inline static void
chi_chru_init( struct unur_gen *gen )
{
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
} /* end of chi_chru_init() */

double
unur_stdgen_sample_chi_chru( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double u,v,z,zz,r;

  if (nu == 1.) {
    while (1) {
      u = uniform();
      v = uniform() * 0.857763884960707;
      z = v / u;
      if (z < 0) continue;
      zz = z * z;
      r = 2.5 - zz;
      if (z < 0.)
	r = r + zz * z / (3. * z);
      if (u < r * 0.3894003915)
	break;
      if (zz > (1.036961043 / u + 1.4))
	continue;
      if (2 * log(u) < (- zz * 0.5 ))
	break;
    }
  }

  else { /* nu > 1 */
    while (1) {
      u = uniform();
      v = uniform() * vd + vm;
      z = v / u;
      if (z < -b)
	continue;
      zz = z * z;
      r = 2.5 - zz;
      if (z < 0.0)
	r = r + zz * z / (3.0 * (z + b));
      if (u < r * 0.3894003915) {
	z += b;
	break;
      }
      if (zz > (1.036961043 / u + 1.4))
	continue;
      if (2. * log(u) < (log(1.0 + z / b) * b * b - zz * 0.5 - z * b)) {
	z += b;
	break;
      }
    }
  }
  /* -X- end of generator code -X- */
  
  return z;

} /* end of unur_stdgen_sample_chi_chru() */

/*---------------------------------------------------------------------------*/
#undef b 
#undef vm
#undef vp
#undef vd
/*---------------------------------------------------------------------------*/
