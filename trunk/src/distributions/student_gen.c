/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      student_gen.c                                                *
 *                                                                           *
 *   Special generators for Student's t distribution                         *
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

inline static void student_trouo_init( struct unur_gen *gen );

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
_unur_stdgen_student_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Student's distribution              */
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

  case 0:  /* Polar Method */     /* Default */
    _unur_cstd_set_sampling_routine( par,gen,unur_stdgen_sample_student_tpol );
    /* no student_tpol_init( gen ) required */
    return 1;

  case 1:  /* Ratio of Uniforms */
    if (par->distr->data.cont.params[0] < 1.) {
      _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
      return 0;
    }
    /* nu >= 1 !!!! */
    _unur_cstd_set_sampling_routine( par,gen,unur_stdgen_sample_student_trouo );
    student_trouo_init( gen );
    return 1;

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default:
    /* no such generator */
    _unur_warning(par->genid,UNUR_ERR_DISTR_GEN,"");
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
 * Student's t Distribution: Polar Method                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from Student's t distribution with  *    
 *               parameters nu > 0.                                          *
 *                                                                           *
 * REFERENCE:  - R.W. Bailey (1994): Polar generation of random variates     *
 *               with the t-distribution,                                    *
 *               Mathematics of Computation 62, 779-781.                     *
 *                                                                           *
 * Implemented by F. Niederl, 1994                                           *
 *****************************************************************************
 *                                                                           *
 * The polar method of Box/Muller for generating Normal variates is adapted  *
 * to the Student-t distribution. The two generated variates are not         *
 * independent and the expected no. of uniforms per variate is 2.5464.       *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*
inline static void student_tpol_init( struct unur_gen *gen )

not required

*/

double
unur_stdgen_sample_student_tpol( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double u,v,w;

  do {
    u = 2. * uniform() - 1.;
    v = 2. * uniform() - 1.;
    w = u * u + v * v;
  } while (w > 1.);

  return(u * sqrt( nu * ( exp(- 2. / nu * log(w)) - 1.) / w));
  /* -X- end of generator code -X- */
} /* end of unur_stdgen_sample_student_tpol() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Student's t Distribution: Ratio of Uniforms                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from Student's t distribution with  *    
 *               parameters nu >= 1.                                         *
 *                                                                           *
 * REFERENCE:  - A.J. Kinderman, J.F. Monahan (1980):                        *
 *               New methods for generating Student's t and gamma variables, *
 *               Computing 25, 369-377.                                      *
 *                                                                           *
 * Implemented by R. Kremer, 1990                                            *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define c       (GEN.gen_param[0])
#define e       (GEN.gen_param[1])
#define p       (GEN.gen_param[2])
#define q       (GEN.gen_param[3])
#define r       (GEN.gen_param[4])
#define vm      (GEN.gen_param[5])
/*---------------------------------------------------------------------------*/

inline static void
student_trouo_init( struct unur_gen *gen )
{
  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = 6;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  if (nu < 1.) {
    _unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
    return;
  }

  r = 1. / nu;
  p = 1. / (1. + r);
  q = -0.25 * (nu + 1.);
  c = 4. * pow(p, q);
  e = 16. / c;
  vm = (nu>1.0) ? sqrt(p+p) * pow( (1.-r)*p, 0.25*(nu-1.) ) : 1.;
  /* -X- end of setup code -X- */
} /* end of student_trouo_init() */

double
unur_stdgen_sample_student_trouo( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double tru,u,v;

  while (1) {

    /* step 1 */
    u = uniform();

    /* step 2 */
    v = uniform();
    v = vm * (v + v - 1.);
    tru = v / u;

    /* step 3 */
    if ( c * u <= 5. - tru * tru) 
      break;
    if (nu >= 3.) 
      if (u * (tru * tru + 3.) >= e) 
	continue;  /* goto 1 */      /* step 4 */
    if ( u <= pow(1. + tru * tru * r, q))
      break;
  }

  return tru;

} /* end of unur_stdgen_sample_student_trouo() */

/*---------------------------------------------------------------------------*/
#undef c
#undef e
#undef p
#undef q
#undef r
#undef vm
/*---------------------------------------------------------------------------*/
