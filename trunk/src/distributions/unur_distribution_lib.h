/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_distribution_lib.h                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros, typedefs and constants and for p.d.f., c.d.f.,    *
 *         etc. of distribtions.                                             *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in distribution source files.                       *
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
#ifndef __UNUR_DISTRIBUTION_LIB_H_SEEN
#define __UNUR_DISTRIBUTION_LIB_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Macros                                                                    */

/* (log of) normalization constant for p.d.f */
#define LOGNORMCONSTANT  params[UNUR_DISTR_MAXPARAMS]
#define NORMCONSTANT     params[UNUR_DISTR_MAXPARAMS]


/* set routine for sampling                                                  */
#if UNUR_DEBUG & UNUR_DB_INFO
#define _unur_cstd_set_sampling_routine(par,gen,routine) \
   do { \
     if ((gen)==NULL) return 1;                       /* test existence only */ \
     (gen)->sample.cont = (routine);                  /* set pointer */ \
     (par)->data.cstd.sample_routine_name = #routine; /* set routine name */ \
   } while (0)
#else
#define _unur_cstd_set_sampling_routine(par,gen,routine) \
   do { \
     if ((gen)==NULL) return 1;                       /* test existence only */ \
     (gen)->sample.cont = (routine);                  /* set pointer */ \
   } while (0)
#endif

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
#endif  /* __UNUR_DISTRIBUTION_LIB_H_SEEN */
/*---------------------------------------------------------------------------*/


