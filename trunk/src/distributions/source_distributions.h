/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_distribution.h                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and prototypes for special generators              *
 *         of distribtions.                                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in distribution source files.                       *
 *                                                                           *
 *                                                                           *
 *   NAMING SCHEME:                                                          *
 *         Includes all distribution source files.                           *
 *         (Not all of these function exist for every distribution!)         *
 *                                                                           *
 *      double _unur_pdf_<distr>(double x, UNUR_DISTR *distr)                *
 *          ... value of PDF at x                                            *
 *      double _unur_dpdf_<distr>(double x, UNUR_DISTR *distr)               *
 *          ... value of derivative of PDF at x                              *
 *      double _unur_cdf_<distr>(double x, UNUR_DISTR *distr)                *
 *          ... value of CDF at x                                            *
 *      double _unur_mode_<distr>(double *params, int n_params)              *
 *          ... mode of distribution                                         *
 *      double _unur_normconstant_<distr>(double *params, int n_params)      *
 *          ... normalization constant of PDF                                *
 *      double _unur_lognormconstant_<distr>(double *params, int n_params)   *
 *          ... log of normalization constant of PDF                         *
 *      int    _unur_stdgen_<distr>_init(                                    *
 *                             UNUR_PAR *parameters, UNUR_GEN *generator)    *
 *          ... initialize new (special) generator for distribution          *
 *      double unur_stdgen_sample_<distr>_<xx>(UNUR_GEN *generator)          *
 *          ... call (special) generator <xx> for distribution               *
 *                                                                           *
 *      double      x          ... argument of PDF                           *
 *      double*     params     ... parameter list for PDF                    *
 *      int         n_params   ... number of parameters (length of array)    *                       
 *                                                                           *
 *      UNUR_DISTR* distr      ... pointer to distribution object            *
 *      UNUR_PAR*   parameters ... pointer to paraters of generator          *
 *      UNUR_GEN*   generator  ... pointer to generator object               *
 *                                                                           *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [1] N.L. Johnson, S. Kotz, and A.W. Kemp                                *
 *       Univariate Discrete Distributions,                                  *
 *       2nd edition                                                         *
 *       John Wiley & Sons, Inc., New York, 1992                             *
 *                                                                           *
 *   [2] N.L. Johnson, S. Kotz, and N. Balakrishnan                          *
 *       Continuous Univariate Distributions,                                *
 *       Volume 1, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1994                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz, and N. Balakrishnan                          *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *   [4] N.L. Johnson, S. Kotz, and N. Balakrishnan                          *
 *       Discrete Multivariate Distributions,                                *
 *       John Wiley & Sons, Inc., New York, 1997                             *
 *                                                                           *
 *   [5] S. Kotz, N. Balakrishnan, and N.L. Johnson                          *
 *       Continuous Multivariate Distributions,                              *
 *       Volume 1: Models and Applications                                   *
 *       John Wiley & Sons, Inc., New York, 2000                             *
 *                                                                           *
 *   [0] S. Kotz and N.L. Johnson                                            *
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
#ifndef __SOURCE_DISTRIBUTIONS_LIB_H_SEEN
#define __SOURCE_DISTRIBUTIONS_LIB_H_SEEN
/*---------------------------------------------------------------------------*/

#include <source_unuran.h>
#include <source_specfunct.h>

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                    Continuous univariate distributions                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/*  Beta distribution  [3; ch.25, p.210]                                     */

/* initialize special generator                                              */
int _unur_stdgen_beta_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_beta_inv( UNUR_GEN *generator );

/* Acceptance/Rejection from log-logistic hats                               */
double _unur_stdgen_sample_beta_bb( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_bc( UNUR_GEN *generator );

/* Stratified Rejection/Patchwork Rejection                                  */
double _unur_stdgen_sample_beta_b00( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_b01( UNUR_GEN *generator );
double _unur_stdgen_sample_beta_b1prs( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Burr family of distributions  [2; ch.12, p.54]                            */

/* initialize special generator                                              */
int _unur_stdgen_burr_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_burr_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Cauchy distribution  [2; ch.16, p.299]                                   */

/* initialize special generator                                              */
int _unur_stdgen_cauchy_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_cauchy_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Chi distribution  [2; ch.18, p.417]                                      */

/* initialize special generator                                              */
int _unur_stdgen_chi_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Ratio of Uniforms with shift                                              */
double _unur_stdgen_sample_chi_chru( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Chisquare distribution  [2; ch.18, p.416]                                */


/*---------------------------------------------------------------------------*/
/* Erlang distribution                                                       */


/*---------------------------------------------------------------------------*/
/*  Exponential distribution  [2; ch.19, p.494]                              */

/* initialize special generator                                              */
int _unur_stdgen_exponential_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_exponential_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Extreme value type I distribution  [3; ch.22, p.2]                       */

/* initialize special generator                                              */
int _unur_stdgen_extremeI_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_extremeI_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Extreme value type II distribution  [3; ch.22, p.2]                      */

/* initialize special generator                                              */
int _unur_stdgen_extremeII_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_extremeII_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Gamma distribution  [2; ch.17, p.337]                                    */

/* initialize special generator                                              */
int _unur_stdgen_gamma_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_gamma_inv( UNUR_GEN *generator );

/* Rejection with log-logistic envelopes                                     */
double _unur_stdgen_sample_gamma_gll( UNUR_GEN *generator );

/* Acceptance Rejection combined with Acceptance Complement */
double _unur_stdgen_sample_gamma_gs( UNUR_GEN *generator );
double _unur_stdgen_sample_gamma_gd( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Generalized inverse Gaussian distribution  [2; ch.15, p.284]              */

/* initialize special generator                                              */
int _unur_stdgen_gig_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Ratio of Uniforms                                                         */
double _unur_stdgen_sample_gig_gigru( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Laplace distribution  [3; ch.24, p.164]                                  */

/* initialize special generator                                              */
int _unur_stdgen_laplace_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_laplace_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Logistic distribution  [3; ch.23, p.115]                                  */

/* initialize special generator                                              */
int _unur_stdgen_logistic_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_logistic_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Lognormal distribution  [2; ch.14, p.208]                                */


/*---------------------------------------------------------------------------*/
/*  Lomax distribution (Pareto distr. of second kind)  [2; ch.20, p.575]     */


/*---------------------------------------------------------------------------*/
/*  Normal distribution  [2; ch.13, p.80]                                    */

/* initialize special generator                                              */
int _unur_stdgen_normal_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_normal_inv( UNUR_GEN *generator );

/* Box-Muller method                                                         */
double _unur_stdgen_sample_normal_bm( UNUR_GEN *generator );

/* Polarmethod with rejection                                                */
double _unur_stdgen_sample_normal_pol( UNUR_GEN *generator );

/* Ratio-of-uniforms method with squeeze                                     */
double _unur_stdgen_sample_normal_quo( UNUR_GEN *generator );

/* "Naive" ratio-of-uniforms method                                          */
double _unur_stdgen_sample_normal_nquo( UNUR_GEN *generator );

/* Ratio-of-uniforms method  with quadratic bounding curves                  */
double _unur_stdgen_sample_normal_leva( UNUR_GEN *generator );

/* Kindermann-Ramage method                                                  */
double _unur_stdgen_sample_normal_kr( UNUR_GEN *generator );

/* Acceptance-complement ratio                                               */
double _unur_stdgen_sample_normal_acr( UNUR_GEN *generator );

/* infamous sum-of-12-uniforms method. NEVER use it!!                        */
double _unur_stdgen_sample_normal_sum( UNUR_GEN *generator );

/*---------------------------------------------------------------------------*/
/*  Pareto distribution (of first kind)  [2; ch.20, p.574]                   */


/*---------------------------------------------------------------------------*/
/* Pearson VI distribution                                                   */


/*---------------------------------------------------------------------------*/
/* Perks distribution                                                        */


/*---------------------------------------------------------------------------*/
/* Planck distribution                                                       */


/*---------------------------------------------------------------------------*/
/*  Power-exponential (Subbotin) distribution  [3; ch.24, p.195]             */

/* initialize special generator                                              */
int _unur_stdgen_powerexponential_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* rejection method for logconcave densities                                 */
double _unur_stdgen_sample_powerexponential_epd( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Rayleigh distribution  [2; ch.18, p.456]                                 */


/*---------------------------------------------------------------------------*/
/* Snedecor's F distribution                                                 */


/*---------------------------------------------------------------------------*/
/* Student's t distribution  [3; ch. 28; p. 362]                             */

/* initialize special generator                                              */
int _unur_stdgen_student_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Polar Method                                                              */
double _unur_stdgen_sample_student_tpol( UNUR_GEN *generator );

/* Ratio of Uniforms                                                         */
double _unur_stdgen_sample_student_trouo( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Slash distribution  [2; ch.12, p.63]                                      */

/* initialize special generator                                              */
int _unur_stdgen_slash_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Ratio of normal and uniform random variates */
double _unur_stdgen_sample_slash_slash( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/*  Triangular distribution  [3; ch.26, p.297]                               */

/* initialize special generator                                              */
int _unur_stdgen_triangular_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_triangular_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Uniform distribution  [3; ch.26, p.276]                                   */

/* initialize special generator                                              */
int _unur_stdgen_uniform_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_uniform_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Weibull distribution  [2; ch.21, p.628]                                   */

/* initialize special generator                                              */
int _unur_stdgen_weibull_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion method                                                          */
double _unur_stdgen_sample_weibull_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                   Continuous multivariate distributions                   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Multinormal distribution  [5; ch.45, p.107]                               */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                     Discrete univariate distributions                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Geometric distribution  [1; ch.5.2, p.201]                                */

/* initialize special generator                                              */
int _unur_stdgen_geometric_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Inversion                                                                 */
int _unur_stdgen_sample_geometric_inv( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Logarithmic distribution  [1; ch.7, p.285]                                */

/* initialize special generator                                              */
int _unur_stdgen_logarithmic_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Acceptance Rejection                                                      */
int _unur_stdgen_sample_logarithmic_lsk( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Negative Binomial distribution  [1; ch.5.1, p.200]                        */

/* initialize special generator                                              */
int _unur_stdgen_negativebinomial_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Compound method                                                           */
int _unur_stdgen_sample_negativebinomial_nbp( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Poisson distribution  [1; ch.4, p.151]                                    */

/* initialize special generator                                              */
int _unur_stdgen_poisson_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Tabulated Inversion combined with Acceptance Complement                   */
int _unur_stdgen_sample_poisson_pdtabl( UNUR_GEN *generator );
int _unur_stdgen_sample_poisson_pdac( UNUR_GEN *generator );

/* Tabulated Inversion combined with Patchwork Rejection                     */
int _unur_stdgen_sample_poisson_pdtabl( UNUR_GEN *generator );
int _unur_stdgen_sample_poisson_pprsc( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/
/* Zipf (or Zeta) distribution  [1; ch.11.20, p.465]                         */

/* initialize special generator                                              */
int _unur_stdgen_zipf_init( UNUR_PAR *parameters, UNUR_GEN *generator );

/* Acceptance Rejection                                                      */
int _unur_stdgen_sample_zipf_zet( UNUR_GEN *generator );


/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                               Macros                                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* set routine for sampling                                                  */
#define _unur_cstd_set_sampling_routine(par,gen,routine) \
   do { \
     if ((gen)==NULL) return 1;                    /* test existence only  */ \
     (gen)->sample.cont = (routine);                 /* set pointer        */ \
     if (par) (par)->data.cstd.sample_routine_name = #routine;  /* set routine name */ \
   } while (0)

#define _unur_dstd_set_sampling_routine(par,gen,routine) \
   do { \
     if ((gen)==NULL) return 1;                    /* test existence only  */ \
     (gen)->sample.discr = (routine);                /* set pointer        */ \
     if (par) (par)->data.dstd.sample_routine_name = #routine;  /* set routine name */ \
   } while (0)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_DISTRIBUTIONS_LIB_H_SEEN */
/*---------------------------------------------------------------------------*/
