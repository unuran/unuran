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
 *      double _unur_pdf_<distr>(double x, double *params, int n_params)     *
 *          ... value of p.d.f. at x                                         *
 *      double _unur_dpdf_<distr>(double x, double *params, int n_params)    *
 *          ... value of derivative of p.d.f. at x                           *
 *      double _unur_cdf_<distr>(double x, double *params, int n_params)     *
 *          ... value of c.d.f. at x                                         *
 *      double _unur_mode_<distr>(double *params, int n_params)              *
 *          ... mode of distribution                                         *
 *      double _unur_normconstant_<distr>(double *params, int n_params)      *
 *          ... normalization constant of p.d.f.                             *
 *      double _unur_lognormconstant_<distr>(double *params, int n_params)   *
 *          ... log of normalization constant of p.d.f.                      *
 *      int    _unur_stdgen_<distr>_init(                                    *
 *                             struct unur_par *par, struct unur_gen *gen)   *
 *          ... initialize new (special) generator for distribution          *
 *      double unur_stdgen_sample_<distr>_<xx>(struct unur_gen *gen)         *
 *          ... call (special) generator <xx> for distribution               *
 *                                                                           *
 *      double  x        ... argument of p.d.f.                              *
 *      double* params   ... parameter list for p.d.f.                       *
 *      int     n_params ... number of parameters (length of parameter array)*
 *                                                                           *
 *      struct unur_par *par ... pointer to paraters of generator            *
 *      struct unur_gen *gen ... pointer to generator object                 *
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

#include <unuran.h>

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*  Beta distribution  [3; ch.25, p.210]                                     */
struct unur_distr *unur_distr_beta(double *params, int n_params);

/* special generators */
int _unur_stdgen_beta_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_beta_bb( struct unur_gen *gen );
double unur_stdgen_sample_beta_bc( struct unur_gen *gen );
/* Acceptance/Rejection from log-logistic hats                               */
double unur_stdgen_sample_beta_b00( struct unur_gen *gen );
double unur_stdgen_sample_beta_b01( struct unur_gen *gen );
double unur_stdgen_sample_beta_b1prs( struct unur_gen *gen );
/* Stratified Rejection/Patchwork Rejection                                  */

/*---------------------------------------------------------------------------*/
/* Burr family of distributions  [2; ch.12, p.54]                            */
struct unur_distr *unur_distr_burr(double *params, int n_params);

/* special generators */
int _unur_stdgen_burr_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_burr_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Cauchy distribution  [2; ch.16, p.299]                                   */
struct unur_distr *unur_distr_cauchy(double *params, int n_params);

/* special generators */
int _unur_stdgen_cauchy_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_cauchy_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Chi distribution  [2; ch.18, p.417]                                      */
struct unur_distr *unur_distr_chi(double *params, int n_params);

/* special generators */
int _unur_stdgen_chi_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_chi_chru( struct unur_gen *gen );
/* Ratio of Uniforms with shift                                              */

/*---------------------------------------------------------------------------*/
/*  Chisquare distribution  [2; ch.18, p.416]                                */
struct unur_distr *unur_distr_chisquare(double *params, int n_params);

double _unur_cdf_chisquare(double x, double *params, int n_params);
/* required for chi^2 tests */

/*---------------------------------------------------------------------------*/
/* Erlang distribution                                                       */
/** TODO **/

/*---------------------------------------------------------------------------*/
/*  Exponential distribution  [2; ch.19, p.494]                              */
struct unur_distr *unur_distr_exponential(double *params, int n_params);

/* special generators */
int _unur_stdgen_exponential_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_exponential_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Gamma distribution  [2; ch.17, p.337]                                    */
struct unur_distr *unur_distr_gamma(double *params, int n_params);

/* special generators */
int _unur_stdgen_gamma_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_gamma_gll( struct unur_gen *gen );
/* Rejection with log-logistic envelopes                                     */
double unur_stdgen_sample_gamma_gs( struct unur_gen *gen );
double unur_stdgen_sample_gamma_gd( struct unur_gen *gen );
/* Acceptance Rejection combined with Acceptance Complement */

/*---------------------------------------------------------------------------*/
/* Generalized inverse Gaussian distribution  [2; ch.15, p.284]              */
struct unur_distr *unur_distr_gig(double *params, int n_params);

/* special generators */
int _unur_stdgen_gig_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_gig_gigru( struct unur_gen *gen );
/* Ratio of Uniforms                                                         */

/*---------------------------------------------------------------------------*/
/*  Laplace distribution  [3; ch.24, p.164]                                  */
struct unur_distr *unur_distr_laplace(double *params, int n_params);

/* special generators */
int _unur_stdgen_laplace_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_laplace_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/* Logistic distribution  [3; ch.23, p.115]                                  */
struct unur_distr *unur_distr_logistic(double *params, int n_params);

/* special generators */
int _unur_stdgen_logistic_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_logistic_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/*  Lognormal distribution  [2; ch.14, p.208]                                */
struct unur_distr *unur_distr_lognormal(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Lomax distribution (Pareto distr. of second kind)  [2; ch.20, p.575]     */
struct unur_distr *unur_distr_lomax(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Normal distribution  [2; ch.13, p.80]                                    */
struct unur_distr *unur_distr_normal( double *params, int n_params );

/* special generators */
 int _unur_stdgen_normal_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_normal_bm( struct unur_gen *gen );
/* Box-Muller method                                                         */
double unur_stdgen_sample_normal_pol( struct unur_gen *gen );
/* Polarmethod with rejection                                                */
double unur_stdgen_sample_normal_quo( struct unur_gen *gen );
/* Ratio-of-uniforms method with squeeze                                     */
double unur_stdgen_sample_normal_nquo( struct unur_gen *gen );
/* "Naive" ratio-of-uniforms method                                          */
double unur_stdgen_sample_normal_leva( struct unur_gen *gen );
/* Ratio-of-uniforms method  with quadratic bounding curves                  */
double unur_stdgen_sample_normal_kr( struct unur_gen *gen );
/* Kindermann-Ramage method                                                  */
double unur_stdgen_sample_normal_acr( struct unur_gen *gen );
/* Acceptance-complement ratio                                               */
double unur_stdgen_sample_normal_sum( struct unur_gen *gen );
/* infamous sum-of-12-uniforms method. NEVER use it!!                        */

/*---------------------------------------------------------------------------*/
/*  Pareto distribution (of first kind)  [2; ch.20, p.574]                   */
struct unur_distr *unur_distr_pareto( double *params, int n_params );

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
struct unur_distr *unur_distr_powerexponential(double *params, int n_params);
/* special generators */
int _unur_stdgen_powerexponential_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_powerexponential_epd( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/*  Rayleigh distribution  [2; ch.18, p.456]                                 */
struct unur_distr *unur_distr_rayleigh(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Snedecor's F distribution                                                 */
/** TODO **/

/*---------------------------------------------------------------------------*/
/* Student's t distribution  [3; ch. 28; p. 362]                             */
struct unur_distr *unur_distr_student(double *params, int n_params);

/* special generators */
int _unur_stdgen_student_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_student_tpol( struct unur_gen *gen );
/* Polar Method                                                              */
double unur_stdgen_sample_student_trouo( struct unur_gen *gen );
/* Ratio of Uniforms                                                         */

/*---------------------------------------------------------------------------*/
/* Slash distribution  [2; ch.12, p.63]                                      */
struct unur_distr *unur_distr_slash(double *params, int n_params);

/* special generators */
int _unur_stdgen_slash_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_slash_slash( struct unur_gen *gen );
/* Ratio of normal and uniform random variates */

/*---------------------------------------------------------------------------*/
/*  Triangular distribution  [3; ch.26, p.297]                               */
struct unur_distr *unur_distr_triangular(double *params, int n_params);

/* special generators */
int _unur_stdgen_triangular_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_triangular_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/* Uniform distribution  [3; ch.26, p.276]                                   */
struct unur_distr *unur_distr_uniform(double *params, int n_params);

/* special generators */
int _unur_stdgen_uniform_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_uniform_inv( struct unur_gen *gen );
/* Inversion method                                                          */


/*---------------------------------------------------------------------------*/
/* Weibull distribution  [2; ch.21, p.628]                                   */
struct unur_distr *unur_distr_weibull(double *params, int n_params);

/* special generators */
int _unur_stdgen_weibull_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_weibull_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_DISTRIBUTIONS_H_SEEN */
/*---------------------------------------------------------------------------*/
