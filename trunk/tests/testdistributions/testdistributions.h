/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: testdistributions.h                                               *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         Make distribution objects for special tests.                      *
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
#ifndef UNURAN_TESTDISTRIBUTIONS_H_SEEN
#define UNURAN_TESTDISTRIBUTIONS_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Special continuous univariate distributions                               */

/*  sawtooth distributions with discontinuous PDF                            */
/*  pdf(x) = |x| - floor(|x|)  (and 1 at integers != 0)                      */
/*  (boundary of domain as parameters)                                       */
UNUR_DISTR *unur_distr_sawtooth_discpdf(const double *params, int n_params);

/*  sawtooth distributions with continuous PDF                               */
/*  pdf(x) = y for y < 0.5 and 1-y otherwise (where y = |x| - floor(|x|))    */ 
/*  (boundary of domain as parameters)                                       */
UNUR_DISTR *unur_distr_sawtooth_contpdf(const double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Continuous univariate distributions without logPDF                        */

/*  Beta distribution  [3; ch.25, p.210]                                     */
UNUR_DISTR *unur_distr_beta_wo_logpdf(const double *params, int n_params);

/* Cauchy distribution  [2; ch.16, p.299]                                    */
UNUR_DISTR *unur_distr_cauchy_wo_logpdf(const double *params, int n_params);

/*  Exponential distribution  [2; ch.19, p.494]                              */
UNUR_DISTR *unur_distr_exponential_wo_logpdf(const double *params, int n_params);

/* Gamma distribution  [2; ch.17, p.337]                                     */
UNUR_DISTR *unur_distr_gamma_wo_logpdf(const double *params, int n_params);

/* Normal distribution  [2; ch.13, p.80]                                     */
UNUR_DISTR *unur_distr_normal_wo_logpdf( const double *params, int n_params );

/*  Power-exponential (Subbotin) distribution  [3; ch.24, p.195]             */
UNUR_DISTR *unur_distr_powerexponential_wo_logpdf(const double *params, int n_params);

/*---------------------------------------------------------------------------*/
/* Specially shaped multivariate distributions                               */

/* Multi-Cauchy distribution where RoU region with r=1 is a ball             */
UNUR_DISTR *unur_distr_multicauchy_RoU_ball( int dim );

/*---------------------------------------------------------------------------*/
/* Multivariate distributions with correlation matrix of AR(1) process       */

/*  Multinormal distribution - AR(1)                                         */
UNUR_DISTR *unur_distr_multinormal_ar1(int dim, const double *mean, double rho);

/*  Multicauchy distribution - AR(1)                                         */
UNUR_DISTR *unur_distr_multicauchy_ar1(int dim, const double *mean, double rho);

/*  Multistudent distribution - AR(1)                                        */
UNUR_DISTR *unur_distr_multistudent_ar1(int dim, double df, const double *mean, double rho);

/*---------------------------------------------------------------------------*/
/* Multivariate distributions with equal off-diagonal entries in             */
/* correlation matrix                                                        */

/*  Multinormal distribution - constant rho                                  */
UNUR_DISTR *unur_distr_multinormal_constantrho(int dim, const double *mean, double rho);

/*  Multicauchy distribution - constant rho                                  */
UNUR_DISTR *unur_distr_multicauchy_constantrho(int dim, const double *mean, double rho);

/*  Multistudent distribution - constant rho                                 */
UNUR_DISTR *unur_distr_multistudent_constantrho(int dim, double df, const double *mean, double rho);

/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_TESTDISTRIBUTIONS_H_SEEN */
/*---------------------------------------------------------------------------*/

