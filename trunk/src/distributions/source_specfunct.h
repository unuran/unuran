/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_specfunct.h                                                *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         prototypes and macros for using special functions like erf(),     *
 *         gamma(), beta(), etc., which are imported from other packages.    *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_distributions.h                           *
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
#ifndef __SOURCE_SPECFUNCT_H_SEEN
#define __SOURCE_SPECFUNCT_H_SEEN
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

/** functions related to beta distribution **/

/* incomplete beta integral */
double incbet(double a, double b, double x);
#define _unur_sf_incomplete_beta(x,a,b)   incbet((a),(b),(x))
#define _HAVE_UNUR_SF_INCOMPLETE_BETA

/** functions related to gamma distribution **/

/* logarithm of gamma function */
double lgam(double x);
#define _unur_sf_ln_gamma(x)   lgam(x)
#define _HAVE_UNUR_SF_LN_GAMMA

/* logarithm of factorial */
#ifdef _HAVE_UNUR_SF_LN_GAMMA
#define _unur_sf_ln_factorial(x)   _unur_sf_ln_gamma((x)+1.)
#define _HAVE_UNUR_SF_LN_FACTORIAL
#endif

/* incomplete gamma function */
double igam(double a, double x);
#define _unur_sf_incomplete_gamma(x,a)  igam((a),(x))
#define _HAVE_UNUR_SF_INCOMPLETE_GAMMA

/** functions related to normal distribution **/

/* cdf of normal distribution */
double ndtr(double x);
#define _unur_sf_cdfnormal(x)   ndtr(x)
#define _HAVE_SF_CDFNORMAL

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_SPECFUNCT_H_SEEN */
/*---------------------------------------------------------------------------*/
