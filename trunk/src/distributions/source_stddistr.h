/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_stddistr.h                                                 *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines identifiers for standard distributions                    *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in ../methods/distr.c and source_distributions.h    *
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
#ifndef __SOURCE_STDDISTR_H_SEEN
#define __SOURCE_STDDISTR_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* indentifiers for standard distributions                                   */

enum {
  UNUR_DISTR_GENERIC  = 0x0u,
  /**                                  pdf   cdf   mode  area  gen   doc     */
  UNUR_DISTR_BETA,                /*    X     X     X     X     X     .      */
  UNUR_DISTR_CAUCHY,              /*    X     X     X     X     X     .      */
  UNUR_DISTR_CHI,                 /*    X     X     X     X     X     .      */
  UNUR_DISTR_CHISQUARE,           /*    X     X     X     X     .     .      */
  UNUR_DISTR_EXPONENTIAL,         /*    X     X     X     X     ?     .      */
  UNUR_DISTR_GAMMA,               /*    X     X     X     X     X     .      */
  UNUR_DISTR_LAPLACE,             /*    X     X     X     X     X     .      */
  UNUR_DISTR_LOGNORMAL,           /*    X     .     .     X     .     .      */
  UNUR_DISTR_LOMAX,               /*    X     X     X     X     .     .      */
  UNUR_DISTR_NORMAL,              /*    X     X     X     X     X     .      */
  UNUR_DISTR_PARETO,              /*    X     X     X     X     .     .      */
  UNUR_DISTR_POWEREXPONENTIAL,    /*    X     .     .     X     .     .      */
  UNUR_DISTR_RAYLEIGH,            /*    X     X     X     X     .     .      */
  UNUR_DISTR_STUDENT,             /*    X     .     X     X     ?     .      */
  UNUR_DISTR_UNIFORM,             /*    X     X     X     X     .     .      */
  UNUR_DISTR_WEIBULL,             /*    X     X     X     X     X     .      */
};

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_STDDISTR_H_SEEN */
/*---------------------------------------------------------------------------*/
