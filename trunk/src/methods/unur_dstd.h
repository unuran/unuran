/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dstd.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DSTD                               *
 *         (wrapper for special generators for                               *
 *         Discrete STanDard distributions)                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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

/*
 METHOD: 
 DSTD (Discrete STandarD distributions)

 DESCRIPTION:
 DSTD is a wrapper for special generator for discrete univariate standard
 distributions. It only works for distributions in the 
 (=>) UNURAN library of standard distributions.

 If any other distribution is provided, or no special generator for the
 given standard distribution is provided, the NULL pointer is returned.
 
 For a distribution more than one special generators (`variants') are possible.
 These are selected by a number. For possible variants see
 (=>) UNURAN library of standard distributions.
 However the following are common to all distributions:
     0                     ... the default generator                      
     UNUR_STDGEN_INVERSION ... the inversion method (if available)         

 Sampling from truncated distributions (which can be constructed by 
 changing the default domain of a distribution by means of an
 (=>) unur_distr_cont_set_domain() call) is not implemented yet.

 Notice that changing the domain of the distribution is not allowed.

 It is possible to change the parameters and the domain of the chosen 
 distribution without building a new generator object.
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

UNUR_PAR *unur_dstd_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for generator. It requires a distribution object 
   for a discrete univariant distribution from the 
   (=>) UNURAN library of standard distributions. 
   Using a truncated distribution is not allowed.
*/

UNUR_GEN *_unur_dstd_init( UNUR_PAR *parameters );
/* Initialize new generator.                                                 */

/** 
    double _unur_dstd_sample( UNUR_GEN *gen );
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

void _unur_dstd_free( UNUR_GEN *generator);
/* Destroy generator object.                                                 */

/*...........................................................................*/

int unur_dstd_set_variant( UNUR_PAR *parameters, unsigned variant );
/* 
   Set variant (special algorithm) for sampling from given distribution.
   For possible variants see (=>) UNURAN library of standard distributions.
   Common variants are `0' for the default generator and
   `UNUR_STDGEN_INVERSION' forthe inversion method (if available).
   If the selected variant number is not implemented, this call has no effect.
*/

int unur_dstd_chg_param( UNUR_GEN *gen, double *params, int n_params );
/*
  Change array of parameters of distribution in given generator object.
  Notice that it is not possible to change the number of parameters.
  This function only copies the given arguments into the array of 
  distribution parameters.
  IMPORTANT: The given parameters are not checked against domain errors;
  in opposition to the (=>) unur_<distr>_new().
*/

/*---------------------------------------------------------------------------*/



