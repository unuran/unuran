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
/*  Beta distribution [3; ch.25, p.210]                                      */
struct unur_distr *unur_distr_beta(double *params, int n_params);

/* special generators */
int _unur_stdgen_beta_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_beta_bb( struct unur_gen *gen );
double unur_stdgen_sample_beta_bc( struct unur_gen *gen );
  /* Acceptance/Rejection from log-logistic hats */

/*---------------------------------------------------------------------------*/
/* Burr ?? distribution                                                      */
/** TODO **/

/*---------------------------------------------------------------------------*/
/*  Cauchy distribution [2; ch.16, p.299]                                    */
struct unur_distr *unur_distr_cauchy(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Chisquare distribution [2; ch.18, p.416]                                 */
struct unur_distr *unur_distr_chisquare(double *params, int n_params);

double _unur_cdf_chisquare(double x, double *params, int n_params);
/* required for chi^2 tests */

/*---------------------------------------------------------------------------*/
/* Erlang distribution                                                       */
/** TODO **/

/*---------------------------------------------------------------------------*/
/*  Exponential distribution [2; ch.19, p.494]                               */
struct unur_distr *unur_distr_exponential(double *params, int n_params);

/* special generators */
int _unur_stdgen_exponential_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_exponential_inv( struct unur_gen *gen );
   /* Inversion method */

/*---------------------------------------------------------------------------*/
/*  Gamma distribution [2; ch.17, p.337]                                     */
struct unur_distr *unur_distr_gamma(double *params, int n_params);

/* special generators */
int _unur_stdgen_gamma_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_gamma_gll( struct unur_gen *gen );
   /* Rejection with log-logistic envelopes */

/*---------------------------------------------------------------------------*/
/* Generalized inverse Gaussian distribution                                 */
/** TODO **/

/*---------------------------------------------------------------------------*/
/*  Laplace distribution [3; ch.24, p.164]                                   */
struct unur_distr *unur_distr_laplace(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Logistic distribution                                                     */
/** TODO **/

/*---------------------------------------------------------------------------*/
/*  Lognormal distribution [2; ch.14, p.208]                                 */
struct unur_distr *unur_distr_lognormal(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Lomax distribution (Pareto distr. of second kind) [2; ch.20, p.575]      */
struct unur_distr *unur_distr_lomax(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Normal distribution [2; ch.13, p.80]                                     */
struct unur_distr *unur_distr_normal( double *params, int n_params );

/* special generators */
 int _unur_stdgen_normal_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_normal_bm( struct unur_gen *gen );
   /* Box-Muller method                                                      */
double unur_stdgen_sample_normal_pol( struct unur_gen *gen );
   /* Polarmethod with rejection                                             */
double unur_stdgen_sample_normal_quo( struct unur_gen *gen );
   /* Ratio-of-uniforms method with squeeze                                  */
double unur_stdgen_sample_normal_nquo( struct unur_gen *gen );
   /* "Naive" ratio-of-uniforms method                                       */
double unur_stdgen_sample_normal_leva( struct unur_gen *gen );
   /* Ratio-of-uniforms method  with quadratic bounding curves               */
double unur_stdgen_sample_normal_kr( struct unur_gen *gen );
   /* Kindermann-Ramage method                                               */
double unur_stdgen_sample_normal_acr( struct unur_gen *gen );
   /* Acceptance-complement ratio                                            */
double unur_stdgen_sample_normal_sum( struct unur_gen *gen );
   /* infamous sum-of-12-uniforms method. NEVER use it!!                     */

/*---------------------------------------------------------------------------*/
/*  Pareto distribution (of first kind) [2; ch.20, p.574]                    */
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
/*  Power-exponential (Subbotin) distribution [3; ch.24, p.195]              */
struct unur_distr *unur_distr_powerexponential(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Rayleigh distribution [2; ch.18, p.456]                                  */
struct unur_distr *unur_distr_rayleigh(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Snedecor's F distribution                                                 */
/** TODO **/

/*---------------------------------------------------------------------------*/
/* Student's t distribution                                                  */
struct unur_distr *unur_distr_student(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Uniform distribution                                                      */
struct unur_distr *unur_distr_uniform(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Weibull distribution                                                      */
/** TODO **/

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_DISTRIBUTIONS_H_SEEN */
/*---------------------------------------------------------------------------*/
