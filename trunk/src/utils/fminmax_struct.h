/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: fminmax_struct.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structure for generic functions                          *
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


/* Structure for generic functions */

struct unur_funct_generic {
  double (*f)(double x, void *params);
  void *params;
};


#define UNUR_FUNCT_GENERIC  double (*) (double, void*)
/*-------------------------------------------------------------*   
 * In order to be able to assign a pointer to a user           *
 * defined function of the form double f(double, double *)     *
 * to the member .f of unur_funct_generic a typecast           *
 * should be performed to avoid compiler warnings :            *
 *                                                             *
 *   double f(double x, double *parameter) {                   *
 *     return (expression_of_x_and_parameter);                 *
 *   }                                                         *
 *                                                             *
 *   main () {                                                 *
 *     struct unur_funct_generic fs;                           *
 *                                                             *
 *     fs.f = (UNUR_FUNCT_GENERIC) f ;                         *
 *   }                                                         *
 *                                                             *
 * In the case of pdf's a similar typecast should              *
 * be performed :                                              *
 *                                                             *
 *     fs.f = (UNUR_FUNCT_GENERIC) DISTR.pdf ;                 *
 *                                                             *
 *-------------------------------------------------------------*/
