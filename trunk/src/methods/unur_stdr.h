/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: stdr.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method STDR                               *
 *         (Simple Transformed Density Rejection with universal bounds)      *
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
   =METHOD  STDR   Simple Transformed Density Rejection

   STDR is based on the transformed density rejection but uses universal 
   inequalities for constructing (universal) hats and squeezes.
   It works for all T-concave distributions with T(x) = -1/sqrt(x).

   It requires the (exact) location of the mode and the area below the given
   p.d.f. The rejection constant is 4 for all T-concave distributions with
   unbounded domain and is less than 4 when the domain is bounded.
   Optionally the c.d.f. at the mode can be given to increase the performance
   of the algorithm. Then the rejection constant is reduced by one half
   and even a universal squeeze can (but need not be) used.
   However using squeezes is not recommended unless the evaluation of the p.d.f.
   is rather expensive.
   
   (The mirror principle is not implemented.)

   If the exact location of the mode is not known, then use an approximation
   and provide the (exact) value of the p.d.f. at the mode by means of the
   unur_stdr_set_pdfatmode() call.

   Instead of the (exact) area below the p.d.f. an upper bound can be used
   (which increases the rejection constant of course). But then the squeeze flag
   must not be set.

   It is even possible to give an upper bound for the p.d.f. and an upper bound
   for the area below the p.d.f. However then the (bound for the) area below the 
   p.d.f. has to be multiplied by the ratio between the upper bound and the
   lower bound of the p.d.f. at the mode. Alternatively the upper bound for the 
   p.d.f. at the mode can be multiplied with the ratio between the upper bound 
   and the lower bound for area below the p.d.f.

   It is possible to change the parameters and the domain of the chosen 
   distribution without building a new generator object using the
   unur_stdr_chg_pdfparams() and unur_stdr_chg_domain() call, respectively.
   But then unur_stdr_chg_pdfarea(), unur_stdr_chg_mode and 
   unur_stdr_chg_cdfatmode() have to used to reset the corresponding figures 
   whenever they have changed.
   If the p.d.f. at the mode has been provided by a 
   unur_stdr_set_pdfatmode() call, additionally unur_stdr_chg_pdfatmode() must 
   be used (otherwise this call is not necessary since then this figure is
   computed directly from the p.d.f.).
   If any of mode, p.d.f. or c.d.f. at the mode, or the area below the mode
   has been changed, then unur_reinit() must be executed.
   (Otherwise the generator produces garbage).
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_stdr_new( UNUR_DISTR *distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_stdr_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
/* Set c.d.f. at mode. 
   When set the performance of the algorithm is increased by factor 2.
   However, when the parameters of the distribution are changed
   (=>) unur_stdr_chg_cdfatmode() has to be used to update this value.
*/

int unur_stdr_set_pdfatmode( UNUR_PAR *parameters, double fmode );
/* Set pdf at mode. if set the p.d.f. at the mode is never changed.          
   This is to avoid additional computations, when the p.d.f. does not
   change when parameters of the distributions vary. 
   It is only useful when the p.d.f. at the mode does not change with
   changing parameters for the distribution.
*/

int unur_stdr_set_verify( UNUR_PAR *parameters, int verify );
/* Turn verifying of algorithm while sampling on/off                         */

int unur_stdr_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
/* Set flag for using universal squeeze (default: off).
   using squeezes is only useful when the evaluation of the p.d.f. is 
   (extremely) expensive.
   Using squeezes is automatically disabled when the c.d.f. at the mode
   is not given (then no universal squeezes exist).
*/

/*...........................................................................*/

int unur_stdr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of distribution in given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of 
   distribution parameters.
   IMPORTANT: The given parameters are not checked against domain errors;
   in opposition to the (=>) unur_<distr>_new() call.
*/

int unur_stdr_chg_domain( UNUR_GEN *generator, double left, double right );
/* Change left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_stdr_chg_mode() is required.
   (There is no domain checking as in the unur_init() call.)
*/

int unur_stdr_chg_mode( UNUR_GEN *generator, double mode );
/* Change mode of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

int unur_stdr_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
/* Change c.d.f. at mode of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

int unur_stdr_chg_pdfatmode( UNUR_GEN *generator, double fmode );
/* Change p.d.f. at mode of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

int unur_stdr_chg_pdfarea( UNUR_GEN *generator, double area );
/* Change area below p.d.f. of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/

