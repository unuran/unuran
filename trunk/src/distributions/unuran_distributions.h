/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_distributions.h                                              *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for PDF, CDF, etc.         *
 *         of distribtions.                                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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

/*---------------------------------------------------------------------------*/
/* Burr family of distributions  [2; ch.12, p.54]                            */
UNUR_DISTR *unur_distr_burr(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Cauchy distribution  [2; ch.16, p.299]                                   */
UNUR_DISTR *unur_distr_cauchy(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Chi distribution  [2; ch.18, p.417]                                      */
UNUR_DISTR *unur_distr_chi(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Chisquare distribution  [2; ch.18, p.416]                                */
UNUR_DISTR *unur_distr_chisquare(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Erlang distribution                                                       */


/*---------------------------------------------------------------------------*/
/*  Exponential distribution  [2; ch.19, p.494]                              */
UNUR_DISTR *unur_distr_exponential(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Extreme value type I distribution  [3; ch.22, p.2]                       */
UNUR_DISTR *unur_distr_extremeI(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Extreme value type II distribution  [3; ch.22, p.2]                      */
UNUR_DISTR *unur_distr_extremeII(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Gamma distribution  [2; ch.17, p.337]                                    */
UNUR_DISTR *unur_distr_gamma(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Generalized inverse Gaussian distribution  [2; ch.15, p.284]              */
UNUR_DISTR *unur_distr_gig(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Laplace distribution  [3; ch.24, p.164]                                  */
UNUR_DISTR *unur_distr_laplace(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Logistic distribution  [3; ch.23, p.115]                                  */
UNUR_DISTR *unur_distr_logistic(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Lognormal distribution  [2; ch.14, p.208]                                */
UNUR_DISTR *unur_distr_lognormal(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Lomax distribution (Pareto distr. of second kind)  [2; ch.20, p.575]     */
UNUR_DISTR *unur_distr_lomax(double *params, int n_params);


/*---------------------------------------------------------------------------*/
/*  Normal distribution  [2; ch.13, p.80]                                    */
UNUR_DISTR *unur_distr_normal( double *params, int n_params );

/*---------------------------------------------------------------------------*/
/*  Pareto distribution (of first kind)  [2; ch.20, p.574]                   */
UNUR_DISTR *unur_distr_pareto( double *params, int n_params );

/*---------------------------------------------------------------------------*/
/* Pearson VI distribution                                                   */


/*---------------------------------------------------------------------------*/
/* Perks distribution                                                        */


/*---------------------------------------------------------------------------*/
/* Planck distribution                                                       */


/*---------------------------------------------------------------------------*/
/*  Power-exponential (Subbotin) distribution  [3; ch.24, p.195]             */
UNUR_DISTR *unur_distr_powerexponential(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Rayleigh distribution  [2; ch.18, p.456]                                 */
UNUR_DISTR *unur_distr_rayleigh(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Snedecor's F distribution                                                 */


/*---------------------------------------------------------------------------*/
/* Student's t distribution  [3; ch. 28; p. 362]                             */
UNUR_DISTR *unur_distr_student(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Slash distribution  [2; ch.12, p.63]                                      */
UNUR_DISTR *unur_distr_slash(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/*  Triangular distribution  [3; ch.26, p.297]                               */
UNUR_DISTR *unur_distr_triangular(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Uniform distribution  [3; ch.26, p.276]                                   */
UNUR_DISTR *unur_distr_uniform(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Weibull distribution  [2; ch.21, p.628]                                   */
UNUR_DISTR *unur_distr_weibull(double *params, int n_params);

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                   Continuous multivariate distributions                   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Multinormal distribution  [5; ch.45, p.107]                               */
UNUR_DISTR *unur_distr_multinormal(int dim, double *mean, double *covar);
/* 
   Creates a distribution object for the multinormal distribution with
   @var{dim} components. @var{mean} is an array of size @var{dim}
   A NULL pointer for @var{mean} is interpreted as the zero
   vector (0,@dots{},0).
   @var{covar} is an array of size @var{dim}x@var{dim} and holds the
   covariance matrix, where the rows of the matrix are stored
   consecutively in this array. The NULL pointer can be used
   instead the identity matrix.

   For standard form of the distribution use NULL for @var{mean} and 
   @var{covar}.
*/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *                     Discrete univariate distributions                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Geometric distribution  [1; ch.5.2, p.201]                                */
UNUR_DISTR *unur_distr_geometric(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Logarithmic distribution  [1; ch.7, p.285]                                */
UNUR_DISTR *unur_distr_logarithmic(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Negative Binomial distribution  [1; ch.5.1, p.200]                        */
UNUR_DISTR *unur_distr_negativebinomial(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Poisson distribution  [1; ch.4, p.151]                                    */
UNUR_DISTR *unur_distr_poisson(double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Zipf (or Zeta) distribution  [1; ch.11.20, p.465]                         */
UNUR_DISTR *unur_distr_zipf(double *params, int n_params);

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_DISTRIBUTIONS_H_SEEN */
/*---------------------------------------------------------------------------*/
