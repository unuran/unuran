/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_student.c                                                  *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  distr: Student distribution or t-distribution [3; ch. 28; p. 362]        *
 *                                                                           *
 *  pdf:       f(x) = ( 1 + (x^2)/nu )^(-(nu+1)/2)                           *
 *  domain:    -infinity < x < infintiy                                      *
 *  constant:  sqrt(nu) * Beta(1/2,nu/2)                                     *
 *             = sqrt(pi*nu) * Gamma(nu/2) / Gamma((nu+1)/2)                 *
 *                                                                           *
 *  parameters:  1                                                           *
 *     0: nu > 0  ... shape                                                  *
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

static const char distr_name[] = "student";

/*---------------------------------------------------------------------------*/
/* parameters */
#define nu  params[0]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the cdf of the distribution ? */
#ifdef HAVE_UNUR_SF_CDFSTUDENT
#  define HAVE_CDF
#else
#  undef  HAVE_CDF
#endif

/* can we compute the area below the pdf ? */
#ifdef HAVE_UNUR_SF_LN_GAMMA
#  define HAVE_AREA
#else
#  undef  HAVE_AREA
#endif

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_student( double x, UNUR_DISTR *distr );
static double _unur_dpdf_student( double x, UNUR_DISTR *distr );
#if HAVE_UNUR_SF_CDFSTUDENT
static double _unur_cdf_student( double x, UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_student( UNUR_DISTR *distr );
#ifdef HAVE_AREA
static int _unur_upd_area_student( UNUR_DISTR *distr );
static double _unur_normconstant_student(double *params, int n_params);
#endif
static int _unur_set_params_student( UNUR_DISTR *distr, double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_student( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  return pow( (1. + x*x/nu), (-nu-1.)*0.5 ) / NORMCONSTANT;
}  /* end of _unur_pdf_student() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_student( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  return ( (-nu-1.)*x/nu * pow( (1. + x*x/nu), (-nu-3.)*0.5 ) / NORMCONSTANT );
} /* end of _unur_dpdf_student() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_student(double x, UNUR_DISTR *distr)
{
  return _unur_sf_cdfstudent(x,DISTR.nu);
} /* end of _unur_cdf_student() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_student( UNUR_DISTR *distr )
{
  DISTR.mode = 0.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return 1;
} /* end of _unur_upd_mode_student() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_AREA

int
_unur_upd_area_student( UNUR_DISTR *distr )
{
  /* normalization constant */
  NORMCONSTANT = _unur_normconstant_student(DISTR.params,DISTR.n_params);
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.area = ( _unur_cdf_student( DISTR.domain[1],distr) 
		 - _unur_cdf_student( DISTR.domain[0],distr) );
  return 1;
#else
  return 0;
#endif

} /* end of _unur_upd_area_student() */

/*---------------------------------------------------------------------------*/

double
_unur_normconstant_student( double *params, int n_params )
{
  return( sqrt(M_PI * nu) * exp(_unur_sf_ln_gamma(0.5*nu) - _unur_sf_ln_gamma(0.5*(nu+1.))) );
} /* end of _unur_normconstant_student() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_set_params_student( UNUR_DISTR *distr, double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return 0; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,0);

  /* check parameter nu */
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0.");
    return 0;
  }

  /* copy parameters for standard form */
  DISTR.nu = nu;

  /* default parameters: none */
  /* copy optional parameters:none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_student() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_student( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_STUDENT;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_student_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_student;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_student; /* pointer to derivative of PDF */
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_student;  /* pointer to CDF               */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_AREA
		 UNUR_DISTR_SET_PDFAREA |
#endif
		 UNUR_DISTR_SET_MODE );

  /* set parameters for distribution */
  if (!_unur_set_params_student(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
#ifdef HAVE_AREA
  NORMCONSTANT = _unur_normconstant_student(DISTR.params,DISTR.n_params);
#else
  NORMCONSTANT = 1.;
#endif

  /* mode and area below p.d.f. */
  DISTR.mode = 0.;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_student;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_student; /* funct for computing mode */
#ifdef HAVE_AREA
  DISTR.upd_area  = _unur_upd_area_student; /* funct for computing area */
#endif

  /* return pointer to object */
  return distr;

} /* end of unur_distr_student() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
