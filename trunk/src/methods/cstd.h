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

   =UP  Methods_for_CONT

   =DESCRIPTION
      CSTD is a wrapper for special generators for continuous
      univariate standard distributions. It only works for
      distributions in the UNURAN library of standard distributions
      (@pxref{Stddist,Standard distributions,Standard distributions}).
      If a distribution object is provided that is build from scratch,
      or no special generator for the given standard distribution is
      provided, the NULL pointer is returned.

      For some distributions more than one special generator
      (@emph{variants}) are possible. These can be choosen by a
      unur_cstd_set_variant() call. For possible variants 
      @pxref{Stddist,Standard distributions,Standard distributions}.
      However the following are common to all distributions:

      @table @code
      @item UNUR_STDGEN_DEFAULT
      the default generator.                      
      @item UNUR_STDGEN_FAST
      the fasted available special generator.
      @item UNUR_STDGEN_INVERSION
      the inversion method (if available).
      @end table
      
      Notice that the variant @code{UNUR_STDGEN_FAST} for a special
      generator might be slower than one of the universal algorithms!
      Additional variants may exist for particular distributions.
      
      Sampling from truncated distributions (which can be constructed by 
      changing the default domain of a distribution by means of
      unur_distr_cont_set_domain() or unur_cstd_chg_truncated() calls)
      is possible but requires the inversion method.
   
      It is possible to change the parameters and the domain of the chosen 
      distribution without building a new generator object.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Common variants for all special generators                                */

#define UNUR_STDGEN_DEFAULT   0        /* default algorithm (don't change 0!)*/
#define UNUR_STDGEN_INVERSION (~0u)    /* inversion method                   */
#define UNUR_STDGEN_FAST      (0)      /* fastest algorithm                  */

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_cstd_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for new generator. It requires a distribution object 
   for a continuous univariant distribution from the 
   UNURAN library of standard distributions 
   (@pxref{Stddist,Standard distributions,Standard distributions}).

   Using a truncated distribution is allowed only if the inversion method
   is available and selected by the unur_cstd_set_variant() call immediately 
   after creating the parameter object.
   Use a unur_distr_cont_set_domain() call to get a truncated distribution.
   To change the domain of a (truncated) distribution use the
   unur_cstd_chg_truncated() call.
*/

/*...........................................................................*/

int unur_cstd_set_variant( UNUR_PAR *parameters, unsigned variant );
/*
   Set variant (special generator) for sampling from a given distribution.
   For possible variants 
   @pxref{Stddist,Standard distributions,Standard distributions}.

   Common variants are @code{UNUR_STDGEN_DEFAULT} for the default generator,
   @code{UNUR_STDGEN_FAST} for (one of the) fastest implemented
   special generators, and @code{UNUR_STDGEN_INVERSION} for the
   inversion method (if available). 
   If the selected variant number is not implemented, this call has no effect.
*/

/*...........................................................................*/

int unur_cstd_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of 
   distribution parameters.

   @emph{IMPORTANT:} The given parameters are not checked against domain errors;
   in opposition to the @command{unur_<distr>_new} calls.
*/

int unur_cstd_chg_truncated( UNUR_GEN *generator, double left, double right );
/* 
   Change left and right border of the domain of the (truncated) distribution.
   This is only possible if the inversion method is used.
   Otherwise this call has no effect.

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call.
*/

/* =END */
/*---------------------------------------------------------------------------*/
