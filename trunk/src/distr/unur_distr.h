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

/*  Beta distribution [3; ch.25, p.210]                                      */
double unur_pdf_beta(double x, double *param, int n_param);
double unur_dpdf_beta(double x, double *param, int n_param);
double unur_cdf_beta(double x, double *param, int n_param);
double unur_mode_beta(double *param, int n_param);
double unur_area_beta(double *param, int n_param);

/* Burr ?? distribution                                                      */
/** TODO **/

/*  Cauchy distribution [2; ch.16, p.299]                                    */
double unur_pdf_cauchy(double x, double *param, int n_param);
double unur_dpdf_cauchy(double x, double *param, int n_param);
double unur_cdf_cauchy(double x, double *param, int n_param);
double unur_mode_cauchy(double *param, int n_param);
double unur_area_cauchy(double *param, int n_param);

/*  Chisquare distribution [2; ch.18, p.416]                                 */
double unur_pdf_chisquare(double x, double *param, int n_param);
double unur_dpdf_chisquare(double x, double *param, int n_param);
double unur_cdf_chisquare(double x, double *param, int n_param);

/* Erlang distribution                                                       */
/** TODO **/

/*  Exponential distribution [2; ch.19, p.494]                               */
double unur_pdf_exponential(double x, double *param, int n_param);
double unur_dpdf_exponential(double x, double *param, int n_param);
double unur_cdf_exponential(double x, double *param, int n_param);
double unur_area_exponential(double *param, int n_param);
double unur_mode_exponential(double *param, int n_param);

/*  Gamma distribution [2; ch.17, p.337]                                     */
double unur_pdf_gamma(double x, double *param, int n_param);
double unur_dpdf_gamma(double x, double *param, int n_param);
double unur_cdf_gamma(double x, double *param, int n_param);
double unur_mode_gamma(double *param, int n_param);
double unur_area_gamma(double *param, int n_param);

/* Generalized inverse Gaussian distribution                                 */
/** TODO **/

/*  Laplace distribution [3; ch.24, p.164]                                   */
double unur_pdf_laplace(double x, double *param, int n_param);
double unur_dpdf_laplace(double x, double *param, int n_param);

/* Logistic distribution                                                     */
/** TODO **/

/*  Lognormal distribution [2; ch.14, p.208]                                 */
double unur_pdf_lognormal(double x, double *param, int n_param);
double unur_dpdf_lognormal(double x, double *param, int n_param);

/*  Lomax distribution (Pareto distr. of second kind) [2; ch.20, p.575]      */
double unur_pdf_lomax(double x, double *param, int n_param);
double unur_dpdf_lomax(double x, double *param, int n_param);

/*  Normal distribution [2; ch.13, p.80]                                     */
double unur_pdf_normal(double x, double *param, int n_param);
double unur_dpdf_normal(double x, double *param, int n_param);
double unur_cdf_normal(double x, double *param, int n_param);
double unur_mode_normal(double *param, int n_param);
double unur_area_normal(double *param, int n_param);

/*  Pareto distribution (of first kind) [2; ch.20, p.574]                    */
double unur_pdf_pareto(double x, double *param, int n_param);
double unur_dpdf_pareto(double x, double *param, int n_param);

/* Pearson VI distribution                                                   */
/** TODO **/

/* Perks distribution                                                        */
/** TODO **/

/* Planck distribution                                                       */
/** TODO **/

/*  Power-exponential (Subbotin) distribution [3; ch.24, p.195]              */
double unur_pdf_powerexponential(double x, double *param, int n_param);
double unur_dpdf_powerexponential(double x, double *param, int n_param);

/*  Rayleigh distribution [2; ch.18, p.456]                                  */
double unur_pdf_rayleigh(double x, double *param, int n_param);
double unur_dpdf_rayleigh(double x, double *param, int n_param);

/* Snedecor's F distribution                                                 */
/** TODO **/

/* Student's t distribution                                                  */
double unur_pdf_student(double x, double *param, int n_param);
double unur_dpdf_student(double x, double *param, int n_param);

/* Uniform distribution                                                      */
double unur_pdf_uniform(double x, double *param, int n_param);
double unur_dpdf_uniform(double x, double *param, int n_param);
double unur_cdf_uniform(double x, double *param, int n_param);
double unur_mode_uniform(double *param, int n_param);
double unur_area_uniform(double *param, int n_param);

/* Weibull distribution                                                      */
/** TODO **/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Prototypes for special functions like erf(), gamma(), beta(), etc.      *
 *   which are imported from other packages.                                 *
 *                                                                           *
 *   We use the package CEPHES/LDOUBLE for computing these functions         *
 *   (available from NETLIB, http://www.netlib.org/cephes/                   *
 *   Copyright 1984 - 1994 by Stephen L. Moshier                             *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* CEPHES library                                                            */

/* cdf of beta(a,b) distribution */
extern long double incbetl(long double a, long double b, long double x);
#define _unur_cdf_beta_ext(x,a,b) ((double)incbetl((long double)(a),(long double)(b),(long double)(x)))

/* cdf of chi^2 distribution with nu degrees of freedom */
extern long double chdtrl(long double df, long double x);
#define _unur_cdf_chisquare_ext(x,nu)  ((double)chdtrl((long double)(nu),(long double)(x)))

/* logarithm of gamma function */
extern long double lgaml(long double x);
#define _unur_gammaln_ext(x)  ((double)(lgaml((long double)(x))))
#define _unur_gammaln(x)      _unur_gammaln_ext(x)

/* cdf of gamma(a,b) distribution */
extern long double gdtrl(long double a, long double b, long double x);
#define _unur_cdf_gamma_ext(x,a,b)    ((double)(gdtrl((long double)(b),(long double)(a),(long double)(x))))

/* cdf of normal distribution */
extern long double ndtrl(long double x);
#define _unur_cdf_normal_ext(x) ((double)(ndtrl((long double)(x))))

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_DISTR_H_SEEN */
/*---------------------------------------------------------------------------*/

