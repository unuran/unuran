/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cdfstudent.c                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         calculates the cdf of the students t distribution                 *
 *                                                                           *
 *****************************************************************************
 
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

#include "unur_specfunct_source.h"

/*---------------------------------------------------------------------------*/

/* cdf for the Student's t distribution */
double _unur_sf_cdfstudent(double x, double nu) {
  double xx;
  if (nu==0) return 0;
  xx=1./(1.+x*x/nu);
  if (x>0) return 1-.5*_unur_sf_incomplete_beta(xx,.5*nu,.5)/_unur_sf_incomplete_beta(1,.5*nu,.5);
  else     return   .5*_unur_sf_incomplete_beta(xx,.5*nu,.5)/_unur_sf_incomplete_beta(1,.5*nu,.5);
}

/*---------------------------------------------------------------------------*/
