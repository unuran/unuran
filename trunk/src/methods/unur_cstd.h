/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cstd.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method CSTD                               *
 *         (wrapper for special generators for                               *
 *         Continuous STanDard distributions)                                *
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
  =METHOD  CSTD   Continuous STandarD distributions
  
  CSTD is a wrapper for special generators for continuous univariate standard
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
  (=>) unur_distr_cont_set_domain() call) is possible but requires the 
  inversion method.
  
  It is possible to change the parameters and the domain of the chosen 
  distribution without building a new generator object.

*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/

UNUR_PAR *unur_cstd_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for new generator. It requires a distribution object 
   for a continuous univariant distribution from the 
   (=>) UNURAN library of standard distributions. 
   Using a truncated distribution is allowed only if the inversion method
   is available and selected by the unur_cstd_set_variant() call immediately 
   after creating the parameter object.
*/

/* Initialize new generator.                                                 */
UNUR_GEN *_unur_cstd_init( UNUR_PAR *parameters );

/** 
    double _unur_cstd_sample( UNUR_GEN *gen );
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

/* Destroy generator object.                                                 */
void _unur_cstd_free( UNUR_GEN *generator);

/*...........................................................................*/

int unur_cstd_set_variant( UNUR_PAR *parameters, unsigned variant );
/* 
   Set variant (special algorithm) for sampling from given distribution.
   For possible variants see (=>) UNURAN library of standard distributions.
   Common variants are `0' for the default generator and
   `UNUR_STDGEN_INVERSION' forthe inversion method (if available).
   If the selected variant number is not implemented, this call has no effect.
*/

int unur_cstd_chg_param( UNUR_GEN *gen, double *params, int n_params );
/*
  Change array of parameters of distribution in given generator object.
  Notice that it is not possible to change the number of parameters.
  This function only copies the given arguments into the array of 
  distribution parameters.
  IMPORTANT: The given parameters are not checked against domain errors;
  in opposition to the (=>) unur_<distr>_new().
*/

int unur_cstd_chg_domain( struct unur_gen *gen, double left, double right );
/* 
   Change left and right border of the domain of the (truncated) distribution.
   This is only possible of the inversion method is used.
   Otherwise this call has no effect.
*/


/* =END */
/*---------------------------------------------------------------------------*/




