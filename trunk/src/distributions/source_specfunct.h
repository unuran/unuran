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
 *   We use the package CEPHES/DOUBLE for computing these functions          *
 *   (available from NETLIB, http://www.netlib.org/cephes/                   *
 *   Copyright 1984 - 1994 by Stephen L. Moshier                             *
 *                                                                           *
 *****************************************************************************/

#ifdef HAVE_LIBMD

/*---------------------------------------------------------------------------*/
/* CEPHES library                                                            */
/*---------------------------------------------------------------------------*/

/** functions related to beta distribution **/

/* incomplete beta integral */
double incbet(double a, double b, double x);
#define _unur_sf_incomplete_beta(x,a,b)   incbet((a),(b),(x))
#define HAVE_UNUR_SF_INCOMPLETE_BETA

/* inverse of incomplete beta integral */
double incbi(double a, double b, double x);
#define _unur_sf_inv_incomplete_beta(x,a,b)  incbi((a),(b),(x))
#define HAVE_UNUR_SF_INV_INCOMPLETE_BETA

/** functions related to gamma distribution **/

/* logarithm of gamma function */
double lgam(double x);
#define _unur_sf_ln_gamma(x)   lgam(x)
#define HAVE_UNUR_SF_LN_GAMMA

/* logarithm of factorial */
#define _unur_sf_ln_factorial(x)   _unur_sf_ln_gamma((x)+1.)
#define HAVE_UNUR_SF_LN_FACTORIAL

/* incomplete gamma function */
double igam(double a, double x);
#define _unur_sf_incomplete_gamma(x,a)  igam((a),(x))
#define HAVE_UNUR_SF_INCOMPLETE_GAMMA

/* inverse of incomplete gamma function */
double igami(double a, double x);
#define _unur_sf_inv_incomplete_gamma(x,a)  igami((a),1.-(x))
#define HAVE_UNUR_SF_INV_INCOMPLETE_GAMMA

/** functions related to normal distribution **/

/* normal distribution function */
double ndtr(double x);
#define _unur_sf_cdfnormal(x)   ndtr(x)
#define HAVE_UNUR_SF_CDFNORMAL

/* inverse of normal distribution function */
double ndtri(double x);
#define _unur_sf_inv_cdfnormal(x)   ndtri(x)
#define HAVE_UNUR_SF_INV_CDFNORMAL

/** functions related to Student's t distribution **/
/* there is no CDF for non-integer degrees of freedom */
#undef HAVE_UNUR_SF_CDFSTUDENT


/*---------------------------------------------------------------------------*/
/* end: CEPHES library                                                       */
/*---------------------------------------------------------------------------*/

#else 

/*---------------------------------------------------------------------------*/
/* Use build in functions (if available)                                     */
/*---------------------------------------------------------------------------*/

/** functions related to beta distribution **/

/* NO incomplete beta integral */
#undef HAVE_UNUR_SF_INCOMPLETE_BETA

/* NO inverse of incomplete beta integral */
#undef HAVE_UNUR_SF_INV_INCOMPLETE_BETA

/** functions related to gamma distribution **/

/* NO logarithm of gamma function */
#undef HAVE_UNUR_SF_LN_GAMMA

/* NO logarithm of factorial */
#undef HAVE_UNUR_SF_LN_FACTORIAL

/* NO incomplete gamma function */
#undef HAVE_UNUR_SF_INCOMPLETE_GAMMA

/* NO inverse of incomplete gamma function */
#undef HAVE_UNUR_SF_INV_INCOMPLETE_GAMMA

/** functions related to normal distribution **/

/* NO normal distribution function */
#undef HAVE_UNUR_SF_CDFNORMAL

/* NO inverse of normal distribution function */
#undef HAVE_UNUR_SF_INV_CDFNORMAL

/** functions related to Student's t distribution **/

/* NO CDF for Stundent's t */
#undef HAVE_UNUR_SF_CDFSTUDENT

#endif

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_SPECFUNCT_H_SEEN */
/*---------------------------------------------------------------------------*/
