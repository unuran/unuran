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
 *         Normalization constants OMITTED!                                  *
 *                                                                           *
 *   NAMING SCHEME:                                                          *
 *         (Not all of these function exist for every distribution!)         *
 *                                                                           *
 *      double unur_pdf_<distr>(double x, double *params, int n_params)      *
 *          ... value of p.d.f. at x                                         *
 *      double unur_dpdf_<distr>(double x, double *params, int n_params)     *
 *          ... value of derivative of p.d.f. at x                           *
 *      double unur_cdf_<distr>(double x, double *params, int n_params)      *
 *          ... value of c.d.f. at x                                         *
 *      double unur_mode_<distr>(double *params, int n_params)               *
 *          ... mode of distribution                                         *
 *      double unur_area_<distr>(double *params, int n_params)               *
 *          ... area below p.d.f. (= normalization constant)                 *
 *                                                                           *
 *      double  x        ... argument of p.d.f.                              *
 *      double* params   ... parameter list for p.d.f.                       *
 *      int     n_params ... number of parameters (length of parameter array)*
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
#ifndef __UNUR_DISTR_H_SEEN
#define __UNUR_DISTR_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>
#include <unur_distr.h>

/*---------------------------------------------------------------------------*/

/** TODO **/
#define NOT_UNIMODAL  0
#define RETURN_NULL   0

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Function prototypes                                                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/*  Beta distribution [3; ch.25, p.210]                                      */
struct unur_distr *unur_distr_beta(double *params, int n_params);
double unur_pdf_beta(double x, double *params, int n_params);
double unur_dpdf_beta(double x, double *params, int n_params);
double unur_cdf_beta(double x, double *params, int n_params);
double unur_mode_beta(double *params, int n_params);
double unur_area_beta(double *params, int n_params);
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
double unur_pdf_cauchy(double x, double *params, int n_params);
double unur_dpdf_cauchy(double x, double *params, int n_params);
double unur_cdf_cauchy(double x, double *params, int n_params);
double unur_mode_cauchy(double *params, int n_params);
double unur_area_cauchy(double *params, int n_params);


/*---------------------------------------------------------------------------*/
/*  Chisquare distribution [2; ch.18, p.416]                                 */
struct unur_distr *unur_distr_chisquare(double *params, int n_params);
double unur_pdf_chisquare(double x, double *params, int n_params);
double unur_dpdf_chisquare(double x, double *params, int n_params);
double unur_cdf_chisquare(double x, double *params, int n_params);


/*---------------------------------------------------------------------------*/
/* Erlang distribution                                                       */
/** TODO **/


/*---------------------------------------------------------------------------*/
/*  Exponential distribution [2; ch.19, p.494]                               */
struct unur_distr *unur_distr_exponential(double *params, int n_params);
double unur_pdf_exponential(double x, double *params, int n_params);
double unur_dpdf_exponential(double x, double *params, int n_params);
double unur_cdf_exponential(double x, double *params, int n_params);
double unur_area_exponential(double *params, int n_params);
double unur_mode_exponential(double *params, int n_params);
/* special generators */
int _unur_stdgen_exponential_init( struct unur_par *par, struct unur_gen *gen );
/* initialize new generator                                                  */
double unur_stdgen_sample_exponential_inv( struct unur_gen *gen );
   /* Inversion method */


/*---------------------------------------------------------------------------*/
/*  Gamma distribution [2; ch.17, p.337]                                     */
struct unur_distr *unur_distr_gamma(double *params, int n_params);
double unur_pdf_gamma(double x, double *params, int n_params);
double unur_dpdf_gamma(double x, double *params, int n_params);
double unur_cdf_gamma(double x, double *params, int n_params);
double unur_mode_gamma(double *params, int n_params);
double unur_area_gamma(double *params, int n_params);
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
double unur_pdf_laplace(double x, double *params, int n_params);
double unur_dpdf_laplace(double x, double *params, int n_params);


/*---------------------------------------------------------------------------*/
/* Logistic distribution                                                     */
/** TODO **/


/*---------------------------------------------------------------------------*/
/*  Lognormal distribution [2; ch.14, p.208]                                 */
struct unur_distr *unur_distr_lognormal(double *params, int n_params);
double unur_pdf_lognormal(double x, double *params, int n_params);
double unur_dpdf_lognormal(double x, double *params, int n_params);


/*---------------------------------------------------------------------------*/
/*  Lomax distribution (Pareto distr. of second kind) [2; ch.20, p.575]      */
struct unur_distr *unur_distr_lomax(double *params, int n_params);
double unur_pdf_lomax(double x, double *params, int n_params);
double unur_dpdf_lomax(double x, double *params, int n_params);


/*---------------------------------------------------------------------------*/
/*  Normal distribution [2; ch.13, p.80]                                     */
struct unur_distr *unur_distr_normal( double *params, int n_params );
double unur_pdf_normal(double x, double *params, int n_params);
double unur_dpdf_normal(double x, double *params, int n_params);
double unur_cdf_normal(double x, double *params, int n_params);
double unur_mode_normal(double *params, int n_params);
double unur_area_normal(double *params, int n_params);
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
double unur_pdf_pareto(double x, double *params, int n_params);
double unur_dpdf_pareto(double x, double *params, int n_params);


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
double unur_pdf_powerexponential(double x, double *params, int n_params);
double unur_dpdf_powerexponential(double x, double *params, int n_params);


/*---------------------------------------------------------------------------*/
/*  Rayleigh distribution [2; ch.18, p.456]                                  */
struct unur_distr *unur_distr_rayleigh(double *params, int n_params);
double unur_pdf_rayleigh(double x, double *params, int n_params);
double unur_dpdf_rayleigh(double x, double *params, int n_params);


/*---------------------------------------------------------------------------*/
/* Snedecor's F distribution                                                 */
/** TODO **/


/*---------------------------------------------------------------------------*/
/* Student's t distribution                                                  */
struct unur_distr *unur_distr_student(double *params, int n_params);
double unur_pdf_student(double x, double *params, int n_params);
double unur_dpdf_student(double x, double *params, int n_params);


/*---------------------------------------------------------------------------*/
/* Uniform distribution                                                      */
struct unur_distr *unur_distr_uniform(double *params, int n_params);
double unur_pdf_uniform(double x, double *params, int n_params);
double unur_dpdf_uniform(double x, double *params, int n_params);
double unur_cdf_uniform(double x, double *params, int n_params);
double unur_mode_uniform(double *params, int n_params);
double unur_area_uniform(double *params, int n_params);


/*---------------------------------------------------------------------------*/
/* Weibull distribution                                                      */
/** TODO **/


/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_DISTR_H_SEEN */
/*---------------------------------------------------------------------------*/


