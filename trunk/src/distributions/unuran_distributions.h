/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_distr.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for p.d.f., c.d.f., etc.   *
 *         of distribtions.                                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         included in all distribution source files.                        *
 *         required for every application of distributions library.          *
 *                                                                           *
 *   COMMENT:                                                                *
 *         Normalization constants may be OMITTED!                           *
 *                                                                           *
 *   NAMING SCHEME:                                                          *
 *         Includes all distribution source files.                           *
 *         (Not all of these function exist for every distribution!)         *
 *                                                                           *
 *      double _unur_pdf_<distr>(double x, UNUR_DISTR *distr)                *
 *          ... value of p.d.f. at x                                         *
 *      double _unur_dpdf_<distr>(double x, UNUR_DISTR *distr)               *
 *          ... value of derivative of p.d.f. at x                           *
 *      double _unur_cdf_<distr>(double x, UNUR_DISTR *distr)                *
 *          ... value of c.d.f. at x                                         *
 *      double _unur_mode_<distr>(double *params, int n_params)              *
 *          ... mode of distribution                                         *
 *      double _unur_normconstant_<distr>(double *params, int n_params)      *
 *          ... normalization constant of p.d.f.                             *
 *      double _unur_lognormconstant_<distr>(double *params, int n_params)   *
 *          ... log of normalization constant of p.d.f.                      *
 *      int    _unur_stdgen_<distr>_init(                                    *
 *                             UNUR_PAR *parameters, UNUR_GEN *generator)    *
 *          ... initialize new (special) generator for distribution          *
 *      double unur_stdgen_sample_<distr>_<xx>(UNUR_GEN *generator)          *
 *          ... call (special) generator <xx> for distribution               *
 *                                                                           *
 *      double      x          ... argument of p.d.f.                        *
 *      double*     params     ... parameter list for p.d.f.                 *
 *      int         n_params   ... number of parameters (length of array)    *                       
 *                                                                           *
 *      UNUR_DISTR* distr      ... pointer to distribution object            *
 *      UNUR_PAR*   parameters ... pointer to paraters of generator          *
 *      UNUR_GEN*   generator  ... pointer to generator object               *
 *                                                                           *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [1] N.L. Johnson, S. Kotz and A.W. Kemp                                 *
 *       Univariate Discrete Distributions,                                  *
 *       2nd edition                                                         *
 *       John Wiley & Sons, Inc., New York, 1992                             *
 *                                                                           *
 *   [2] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 1, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1994                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *   [4] S. Kotz and N.L. Johnson                                            *
 *       Encyclopedia of Statistical Sciences                                *
 *       Volumes 1-9                                                         *
 *       John Wiley & Sons, Inc., New York, 1982-1988                        *
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
#ifndef __UNURAN_DISTRIBUTIONS_H_SEEN
#define __UNURAN_DISTRIBUTIONS_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_stddistr.h>

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                    Continuous univariate distributions                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/*  Beta distribution  [3; ch.25, p.210]                                     */
UNUR_DISTR *unur_distr_beta(double *params, int n_params);

/* special generators */
int _unur_stdgen_beta_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_beta_bb( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_bc( UNUR_GEN *generator );
/* Acceptance/Rejection from log-logistic hats                               */
double _unur_stdgen_sample_beta_b00( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_b01( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_b1prs( UNUR_GEN *generator );
/* Stratified Rejection/Patchwork Rejection                                  */

/*---------------------------------------------------------------------------*/
/* Burr family of distributions  [2; ch.12, p.54]                            */
UNUR_DISTR *unur_distr_burr(double *params, int n_params);

/* special generators */
int _unur_stdgen_burr_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_burr_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Cauchy distribution  [2; ch.16, p.299]                                   */
UNUR_DISTR *unur_distr_cauchy(double *params, int n_params);

/* special generators */
int _unur_stdgen_cauchy_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_cauchy_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Chi distribution  [2; ch.18, p.417]                                      */
UNUR_DISTR *unur_distr_chi(double *params, int n_params);

/* special generators */
int _unur_stdgen_chi_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_chi_chru( UNUR_GEN *generator );
/* Ratio of Uniforms with shift                                              */

/*---------------------------------------------------------------------------*/
/*  Chisquare distribution  [2; ch.18, p.416]                                */
UNUR_DISTR *unur_distr_chisquare(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Erlang distribution                                                       */
/** TODO **/

/*---------------------------------------------------------------------------*/
/*  Exponential distribution  [2; ch.19, p.494]                              */
UNUR_DISTR *unur_distr_exponential(double *params, int n_params);

/* special generators */
int _unur_stdgen_exponential_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_exponential_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Extreme value type I distribution  [3; ch.22, p.2]                       */
UNUR_DISTR *unur_distr_extremeI(double *params, int n_params);

/* special generators */
int _unur_stdgen_extremeI_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_extremeI_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Extreme value type II distribution  [3; ch.22, p.2]                      */
UNUR_DISTR *unur_distr_extremeII(double *params, int n_params);

/* special generators */
int _unur_stdgen_extremeII_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_extremeII_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Gamma distribution  [2; ch.17, p.337]                                    */
UNUR_DISTR *unur_distr_gamma(double *params, int n_params);

/* special generators */
int _unur_stdgen_gamma_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_gamma_gll( UNUR_GEN *generator );
/* Rejection with log-logistic envelopes                                     */
double _unur_stdgen_sample_gamma_gs( UNUR_GEN *generator );
double _unur_stdgen_sample_gamma_gd( UNUR_GEN *generator );
/* Acceptance Rejection combined with Acceptance Complement */

/*---------------------------------------------------------------------------*/
/* Generalized inverse Gaussian distribution  [2; ch.15, p.284]              */
UNUR_DISTR *unur_distr_gig(double *params, int n_params);

/* special generators */
int _unur_stdgen_gig_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_gig_gigru( UNUR_GEN *generator );
/* Ratio of Uniforms                                                         */

/*---------------------------------------------------------------------------*/
/*  Laplace distribution  [3; ch.24, p.164]                                  */
UNUR_DISTR *unur_distr_laplace(double *params, int n_params);

/* special generators */
int _unur_stdgen_laplace_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_laplace_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/* Logistic distribution  [3; ch.23, p.115]                                  */
UNUR_DISTR *unur_distr_logistic(double *params, int n_params);

/* special generators */
int _unur_stdgen_logistic_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_logistic_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Lognormal distribution  [2; ch.14, p.208]                                */
UNUR_DISTR *unur_distr_lognormal(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Lomax distribution (Pareto distr. of second kind)  [2; ch.20, p.575]     */
UNUR_DISTR *unur_distr_lomax(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Normal distribution  [2; ch.13, p.80]                                    */
UNUR_DISTR *unur_distr_normal( double *params, int n_params );

/* special generators */
 int _unur_stdgen_normal_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_normal_bm( UNUR_GEN *generator );
/* Box-Muller method                                                         */
double _unur_stdgen_sample_normal_pol( UNUR_GEN *generator );
/* Polarmethod with rejection                                                */
double _unur_stdgen_sample_normal_quo( UNUR_GEN *generator );
/* Ratio-of-uniforms method with squeeze                                     */
double _unur_stdgen_sample_normal_nquo( UNUR_GEN *generator );
/* "Naive" ratio-of-uniforms method                                          */
double _unur_stdgen_sample_normal_leva( UNUR_GEN *generator );
/* Ratio-of-uniforms method  with quadratic bounding curves                  */
double _unur_stdgen_sample_normal_kr( UNUR_GEN *generator );
/* Kindermann-Ramage method                                                  */
double _unur_stdgen_sample_normal_acr( UNUR_GEN *generator );
/* Acceptance-complement ratio                                               */
double _unur_stdgen_sample_normal_sum( UNUR_GEN *generator );
/* infamous sum-of-12-uniforms method. NEVER use it!!                        */

/*---------------------------------------------------------------------------*/
/*  Pareto distribution (of first kind)  [2; ch.20, p.574]                   */
UNUR_DISTR *unur_distr_pareto( double *params, int n_params );

/*---------------------------------------------------------------------------*/
/* Pearson VI distribution                                                   */
/** TODO **/

/*---------------------------------------------------------------------------*/
/* Perks distribution                                                        */
/** TODO **/

/*---------------------------------------------------------------------------*/
/* Planck distribution                                                       */
/** TODO **/

/*---------------------------------------------------------------------------*/
/*  Power-exponential (Subbotin) distribution  [3; ch.24, p.195]             */
UNUR_DISTR *unur_distr_powerexponential(double *params, int n_params);
/* special generators */
int _unur_stdgen_powerexponential_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_powerexponential_epd( UNUR_GEN *generator );

/*---------------------------------------------------------------------------*/
/*  Rayleigh distribution  [2; ch.18, p.456]                                 */
UNUR_DISTR *unur_distr_rayleigh(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Snedecor's F distribution                                                 */
/** TODO **/

/*---------------------------------------------------------------------------*/
/* Student's t distribution  [3; ch. 28; p. 362]                             */
UNUR_DISTR *unur_distr_student(double *params, int n_params);

/* special generators */
int _unur_stdgen_student_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_student_tpol( UNUR_GEN *generator );
/* Polar Method                                                              */
double _unur_stdgen_sample_student_trouo( UNUR_GEN *generator );
/* Ratio of Uniforms                                                         */

/*---------------------------------------------------------------------------*/
/* Slash distribution  [2; ch.12, p.63]                                      */
UNUR_DISTR *unur_distr_slash(double *params, int n_params);

/* special generators */
int _unur_stdgen_slash_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_slash_slash( UNUR_GEN *generator );
/* Ratio of normal and uniform random variates */

/*---------------------------------------------------------------------------*/
/*  Triangular distribution  [3; ch.26, p.297]                               */
UNUR_DISTR *unur_distr_triangular(double *params, int n_params);

/* special generators */
int _unur_stdgen_triangular_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_triangular_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/* Uniform distribution  [3; ch.26, p.276]                                   */
UNUR_DISTR *unur_distr_uniform(double *params, int n_params);

/* special generators */
int _unur_stdgen_uniform_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_uniform_inv( UNUR_GEN *generator );
/* Inversion method                                                          */


/*---------------------------------------------------------------------------*/
/* Weibull distribution  [2; ch.21, p.628]                                   */
UNUR_DISTR *unur_distr_weibull(double *params, int n_params);

/* special generators */
int _unur_stdgen_weibull_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
double _unur_stdgen_sample_weibull_inv( UNUR_GEN *generator );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                     Discrete univariate distributions                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Geometric distribution  [1; ch.5.2, p.201]                                */
UNUR_DISTR *unur_distr_geometric(double *params, int n_params);

/* special generators */
int _unur_stdgen_geometric_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
int _unur_stdgen_sample_geometric_inv( UNUR_GEN *generator );
/* Inversion                                                                 */

/*---------------------------------------------------------------------------*/
/* Logarithmic distribution  [1; ch.7, p.285]                                */
UNUR_DISTR *unur_distr_logarithmic(double *params, int n_params);

/* special generators */
int _unur_stdgen_logarithmic_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
int _unur_stdgen_sample_logarithmic_lsk( UNUR_GEN *generator );
/* Acceptance Rejection                                                      */

/*---------------------------------------------------------------------------*/
/* Negative Binomial distribution  [1; ch.5.1, p.200]                        */
UNUR_DISTR *unur_distr_negativebinomial(double *params, int n_params);

/* special generators */
int _unur_stdgen_negativebinomial_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
int _unur_stdgen_sample_negativebinomial_nbp( UNUR_GEN *generator );
/* Compound method                                                           */

/*---------------------------------------------------------------------------*/
/* Poisson distribution  [1; ch.4, p.151]                                    */
UNUR_DISTR *unur_distr_poisson(double *params, int n_params);

/* special generators */
int _unur_stdgen_poisson_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
int _unur_stdgen_sample_poisson_pdtabl( UNUR_GEN *generator );
int _unur_stdgen_sample_poisson_pdac( UNUR_GEN *generator );
/* Tabulated Inversion combined with Acceptance Complement                   */
int _unur_stdgen_sample_poisson_pdtabl( UNUR_GEN *generator );
int _unur_stdgen_sample_poisson_pprsc( UNUR_GEN *generator );
/* Tabulated Inversion combined with Patchwork Rejection                     */

/*---------------------------------------------------------------------------*/
/* Zipf (or Zeta) distribution  [1; ch.11.20, p.465]                         */
UNUR_DISTR *unur_distr_zipf(double *params, int n_params);

/* special generators */
int _unur_stdgen_zipf_init( UNUR_PAR *parameters, UNUR_GEN *generator );
/* initialize new generator                                                  */
int _unur_stdgen_sample_zipf_zet( UNUR_GEN *generator );
/* Acceptance Rejection                                                      */

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_DISTRIBUTIONS_H_SEEN */
/*---------------------------------------------------------------------------*/
