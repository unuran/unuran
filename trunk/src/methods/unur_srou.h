/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: srou.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method SROU                               *
 *         (Simple universal generator, ratio-of-uniforms method)            *
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
   =METHOD  SROU   Simple Ratio-Of-Uniforms method

   SROU is based on the ratio-of-uniforms method but uses universal 
   inequalities for constructing a (universal) bounding rectangle.
   It works for all T-concave distributions with T(x) = -1/sqrt(x).

   It requires the p.d.f., the (exact) location of the mode and the area below 
   the given p.d.f. The rejection constant is 4 for all T-concave distributions.
   Optionally the c.d.f. at the mode can be given to increase the performance
   of the algorithm. Then the rejection constant is reduced to 2 and even a
   universal squeeze can (but need not be) used.
   A way to increase the performence of the algorithm when the c.d.f. at the
   mode is not provided is the usae of the mirror principle.
   However using squeezes and using the mirror principle is not recommended
   in general (see below).

   If the exact location of the mode is not known, then use the approximate
   location and provide the (exact) value of the p.d.f. at the mode by means
   of the unur_stdr_set_pdfatmode() call.

   If the (exact) area below the p.d.f. is not known, then an upper bound can be
   used instead (which of course increases the rejection constant). 
   But then the squeeze flag must not be set.

   It is even possible to give an upper bound for the p.d.f. only.
   However then the (upper bound for the) area below the p.d.f. has to be 
   multiplied by the ratio between the upper bound and the lower bound of the 
   p.d.f. at the mode.
   Again seting the squeeze flag is not allowed.

   It is possible to change the parameters and the domain of the chosen 
   distribution without building a new generator object using the
   unur_srou_chg_pdfparams() and unur_srou_chg_domain() call, respectively.
   But then unur_srou_chg_pdfarea(), unur_srou_chg_mode and 
   unur_srou_chg_cdfatmode() have to used to reset the corresponding figures 
   whenever they have changed.
   If the p.d.f. at the mode has been provided by a 
   unur_srou_set_pdfatmode() call, additionally unur_srou_chg_pdfatmode() must 
   be used (otherwise this call is not necessary since then this figure is
   computed directly from the p.d.f.).
   If any of mode, p.d.f. or c.d.f. at the mode, or the area below the mode
   has been changed, then unur_reinit() must be executed.
   (Otherwise the generator produces garbage).
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_srou_new( UNUR_DISTR *distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_srou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
/* Set c.d.f. at mode. 
   When set the performance of the algorithm is increased by factor 2.
   However, when the parameters of the distribution are changed
   (=>) unur_srou_chg_cdfatmode() has to be used to update this value.
*/

int unur_srou_set_pdfatmode( UNUR_PAR *parameters, double fmode );
/* Set pdf at mode. if set the p.d.f. at the mode is never changed.          
   This is to avoid additional computations, when the p.d.f. does not
   change when parameters of the distributions vary. 
   It is only useful when the p.d.f. at the mode does not change with
   changing parameters for the distribution.
*/

int unur_srou_set_verify( UNUR_PAR *parameters, int verify );
/* Turn verifying of algorithm while sampling on/off                         */

int unur_srou_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
/* Set flag for using universal squeeze (default: off).
   using squeezes is only useful when the evaluation of the p.d.f. is 
   (extremely) expensive.
   Using squeezes is automatically disabled when the c.d.f. at the mode
   is not given (then no universal squeezes exist).
*/

int unur_srou_set_usemirror( UNUR_PAR *parameters, int usemirror );
/* Set flag for using mirror principle (default: off)                        
   Using the mirror principle is only useful when the c.d.f. at the
   mode is not known and the evaluation of the p.d.f. is rather cheap compared
   to the marginal generation time of the underlying uniform random
   number generator.
   It is automatically disabled when the c.d.f. at the mode is given.
   (Then there is no necessity to use the mirror principle. However disabling
   is only done during the initialization step but not at a re-initialization
   step.)
*/

/*...........................................................................*/

int unur_srou_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of distribution in given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of 
   distribution parameters.
   IMPORTANT: The given parameters are not checked against domain errors;
   in opposition to the (=>) unur_<distr>_new() call.
*/

int unur_srou_chg_domain( UNUR_GEN *generator, double left, double right );
/* Change left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_srou_chg_mode() is required.
   (There is no domain checking as in the unur_init() call.)
*/

int unur_srou_chg_mode( UNUR_GEN *generator, double mode );
/* Change mode of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

int unur_srou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
/* Change c.d.f. at mode of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

int unur_srou_chg_pdfatmode( UNUR_GEN *generator, double fmode );
/* Change p.d.f. at mode of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

int unur_srou_chg_pdfarea( UNUR_GEN *generator, double area );
/* Change area below p.d.f. of distribution.
   unur_reinit() must be executed before sampling from the generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/


