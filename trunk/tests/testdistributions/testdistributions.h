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

/*  Multinormal distribution (corr-matrix from AR1 process)                  */
UNUR_DISTR *unur_distr_multinormal_ar1(int dim, const double *mean, double rho);

/*  Multicauchy distribution (corr-matrix from AR1 process)                  */
UNUR_DISTR *unur_distr_multicauchy_ar1(int dim, const double *mean, double rho);

/*  Multistudent distribution (corr-matrix from AR1 process)                  */
UNUR_DISTR *unur_distr_multistudent_ar1(int dim, double df, const double *mean, double rho);

/*  Multinormal distribution (corr-matrix with equal off-diagonal elements)  */
UNUR_DISTR *unur_distr_multinormal_constant_rho(int dim, const double *mean, double rho);

/*  Multicauchy distribution (corr-matrix with equal off-diagonal elements)  */
UNUR_DISTR *unur_distr_multicauchy_constant_rho(int dim, const double *mean, double rho);

/*  Multistudent distribution (corr-matrix with equal off-diagonal elements)  */
UNUR_DISTR *unur_distr_multistudent_constant_rho(int dim, double df, const double *mean, double rho);

/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_TESTDISTRIBUTIONS_H_SEEN */
/*---------------------------------------------------------------------------*/

