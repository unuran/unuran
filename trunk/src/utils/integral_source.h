/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: integral_source.h                                                 *
 *                                                                           *
 *   Routines for integral computations.                                     *
 *                                                                           *
 *****************************************************************************
 
 *****************************************************************************
 *                                                                           *
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

int
_unur_integral_evaluate_single (int dim, struct unur_funct_vgeneric *f, UNUR_GEN *gen, 
                                long n_samples, double *integral);
/* Integrating f(x) in R^dim with x being distributed according to generator gen */

int
_unur_integral_evaluate (int dim, struct unur_funct_vgeneric *f, UNUR_GEN *gen, 
			 long n_samples, long n_repetitions, double *exact_value,
			 double *mean, double *variance, double *rmse);
/* Integrating f(x) in R^dim with x being distributed according to generator gen */

