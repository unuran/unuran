/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_stddistr.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines identifiers for standard distributions                    *
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
#ifndef UNUR_STDDISTR_H_SEEN
#define UNUR_STDDISTR_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* indentifiers for standard distributions                                   */

enum {
  UNUR_DISTR_GENERIC  = 0x0u,

  UNUR_DISTR_CORDER,              /* order statistics                        */

  /**                                  pdf   cdf   mode  area  gen   doc     */
  UNUR_DISTR_BETA,                /*    X     X     X     X     X     .      */
  UNUR_DISTR_CAUCHY,              /*    X     X     X     X     X     .      */
  UNUR_DISTR_CHI,                 /*    X     X     X     X     X     .      */
  UNUR_DISTR_CHISQUARE,           /*    X     X     X     X     .     .      */
  UNUR_DISTR_EPANECHNIKOV,        /*    .     .     .     .     .     .      */
  UNUR_DISTR_EXPONENTIAL,         /*    X     X     X     X     X     .      */
  UNUR_DISTR_EXTREME_I,           /*    X     X     X     X     X     .      */
  UNUR_DISTR_EXTREME_II,          /*    X     X     X     X     X     .      */
  UNUR_DISTR_GAMMA,               /*    X     X     X     X     X     .      */
  UNUR_DISTR_GIG,                 /*    X     X     .     .     X     .      */
  UNUR_DISTR_LAPLACE,             /*    X     X     X     X     X     .      */
  UNUR_DISTR_LOGISTIC,            /*    X     X     X     X     X     .      */
  UNUR_DISTR_LOGNORMAL,           /*    X     .     .     X     .     .      */
  UNUR_DISTR_LOMAX,               /*    X     X     X     X     .     .      */
  UNUR_DISTR_NORMAL   = 0x0100u,  /*    X     X     X     X     X     .      */
   UNUR_DISTR_GAUSSIAN = 0x0100u, /*    same as NORMAL                       */
  UNUR_DISTR_PARETO,              /*    X     X     X     X     .     .      */
  UNUR_DISTR_POWEREXPONENTIAL,    /*    X     X     X     X     X     .      */
  UNUR_DISTR_RAYLEIGH,            /*    X     X     X     X     .     .      */
  UNUR_DISTR_SLASH,               /*    X     .     X     X     X     .      */
  UNUR_DISTR_STUDENT,             /*    X    (X)    X     X     X     .      */
  UNUR_DISTR_TRIANGULAR,          /*    X     X     X     X     X     .      */
  UNUR_DISTR_UNIFORM = 0x0200u,   /*    X     X     X     X     .     .      */
   UNUR_DISTR_BOXCAR = 0x0200u,   /*    same as UNIFORM                      */
  UNUR_DISTR_WEIBULL,             /*    X     X     X     X     X     .      */

  UNUR_DISTR_BURR_I,              /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_II,             /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_III,            /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_IV,             /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_V,              /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_VI,             /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_VII,            /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_VIII,           /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_IX,             /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_X,              /*    .     X     .     .     X     .      */
  UNUR_DISTR_BURR_XI,             /*    .     X     .     .     .     .      */
  UNUR_DISTR_BURR_XII,            /*    .     X     .     .     X     .      */

  /**                                  pmf   cdf   mode  sum   gen   doc     */
  UNUR_DISTR_BINOMIAL,            /*    .     .     .     .     .     .      */
  UNUR_DISTR_GEOMETRIC,           /*    X     X     X     X     X     .      */
  UNUR_DISTR_HYPERGEOMETRIC,      /*    .     .     .     .     .     .      */
  UNUR_DISTR_LOGARITHMIC,         /*    X     .     X     X     X     .      */
  UNUR_DISTR_NEGATIVEBINOMIAL,    /*    X     .     X     X     X     .      */
  UNUR_DISTR_POISSON,             /*    X     X     X     X     X     .      */
  UNUR_DISTR_ZIPF,                /*    X     .     X     .     X     .      */

  /**                                  pdf                           doc     */
  UNUR_DISTR_MNORMAL,             /*    X                             .      */
  UNUR_DISTR_COPULA,              /*    .                             .      */

  /**                                                                        */
  UNUR_DISTR_MCORRELATION         /*                                         */

};

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_STDDISTR_H_SEEN */
/*---------------------------------------------------------------------------*/
