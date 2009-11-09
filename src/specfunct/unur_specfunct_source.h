/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_specfunct_source.h                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         prototypes and macros for using special functions like erf(),     *
 *         gamma(), beta(), etc., which are imported from other packages.    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
#ifndef UNUR_SPECFUNCT_SOURCE_H_SEEN
#define UNUR_SPECFUNCT_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Prototypes for special functions like erf(), gamma(), beta(), etc.      *
 *   which are imported from other packages.                                 *
 *                                                                           *
 *   We use the package CEPHES/DOUBLE for computing these functions          *
 *   (available from NETLIB, http://www.netlib.org/cephes/                   *
 *   Copyright 1984 - 1994 by Stephen L. Moshier                             *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Routines from the CEPHES library.                                         */
/*---------------------------------------------------------------------------*/

/** functions related to beta distribution **/

/* incomplete beta integral */
double _unur_cephes_incbet(double a, double b, double x);
#define _unur_sf_incomplete_beta(x,a,b)   _unur_cephes_incbet((a),(b),(x))

/** functions related to gamma distribution **/

/* logarithm of gamma function */
double _unur_cephes_lgam(double x);
#define _unur_sf_ln_gamma(x)   _unur_cephes_lgam(x)

/* logarithm of factorial */
#define _unur_sf_ln_factorial(x)   _unur_sf_ln_gamma((x)+1.)

/* incomplete gamma function */
double _unur_cephes_igam(double a, double x);
#define _unur_sf_incomplete_gamma(x,a)  _unur_cephes_igam((a),(x))

/** functions related to normal distribution **/

/* normal distribution function */
double _unur_cephes_ndtr(double x);
#define _unur_sf_cdfnormal(x)   _unur_cephes_ndtr(x)

/* inverse of normal distribution function */
double _unur_cephes_ndtri(double x);
#define _unur_sf_inv_cdfnormal(x)   _unur_cephes_ndtri(x)

/*---------------------------------------------------------------------------*/
/* end: CEPHES library                                                       */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Routines from the Rmath library (R project).                              */
/*---------------------------------------------------------------------------*/

#ifdef HAVE_LIBRMATH

#define HAVE_BESSEL_K 1

/* modified Bessel function K_nu of second kind (AKA third kind)             */
double _unur_sf_bessel_k(double x, double nu);

/* rescaled modified Bessel function K_nu of second kind (AKA third kind):   */
/*   exp(x) * bessel_k(x,nu)                                                 */
double _unur_sf_bessel_k_expo(double x, double nu);
#endif

/*---------------------------------------------------------------------------*/
/* end: Rmath library (R project)                                            */
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Replacement for missing (system) functions                              *
 *                                                                           *
 *****************************************************************************/

#if !HAVE_DECL_LOG1P
/* log(1+x) */
/* (replacement for missing C99 function log1p) */
double _unur_log1p(double x);
#endif

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_SPECFUNCT_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
