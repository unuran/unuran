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
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
   =METHOD  DSTD    Discrete STandarD distributions
  
   =UP  Methods_for_DISCR

   =REQUIRED standard distribution from UNURAN library
      (@pxref{Stddist,,Standard distributions}).

   =SPEED Set-up: fast, Sampling: depends on distribution and generator

   =DESCRIPTION
      DSTD is a wrapper for special generators for discrete univariate
      standard distributions. It only works for distributions in the 
      UNURAN library of standard distributions
      (@pxref{Stddist,,Standard distributions}).
      If a distribution object is provided that is build from scratch,
      or no special generator for the given standard distribution is
      provided, the NULL pointer is returned.

      For some distributions more than one special generator
      is possible. 
      
   =HOWTOUSE
      Create a distribution object for a standard distribution
      from the UNURAN library (@pxref{Stddist,,Standard distributions}).
      For some distributions more than one special generator
      (@emph{variants}) is possible. These can be choosen by a
      unur_dstd_set_variant() call. For possible variants 
      @xref{Stddist,,Standard distributions}.
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
      unur_distr_discr_set_domain() call)
      is possible but requires the inversion method.

      It is possible to change the parameters and the domain of the chosen 
      distribution without building a new generator object
      by means of unur_dstd_chg_pmfparams().

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/

UNUR_PAR *unur_dstd_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for new generator. It requires a distribution object 
   for a discrete univariant distribution from the 
   UNURAN library of standard distributions 
   (@pxref{Stddist,,Standard distributions}).

   Using a truncated distribution is allowed only if the inversion method
   is available and selected by the unur_dstd_set_variant() call immediately 
   after creating the parameter object.
   Use a unur_distr_discr_set_domain() call to get a truncated
   distribution.
*/

/*...........................................................................*/

int unur_dstd_set_variant( UNUR_PAR *parameters, unsigned variant );
/* 
   Set variant (special generator) for sampling from a given distribution.
   For possible variants 
   @pxref{Stddist,,Standard distributions}.

   Common variants are @code{UNUR_STDGEN_DEFAULT} for the default generator,
   @code{UNUR_STDGEN_FAST} for (one of the) fastest implemented
   special generators, and @code{UNUR_STDGEN_INVERSION} for the
   inversion method (if available). 
   If the selected variant number is not implemented, this call has no effect.
*/

int unur_dstd_chg_pmfparams( UNUR_GEN *gen, double *params, int n_params );
/*
   Change array of parameters of the distribution in a given generator
   object. If the given parameters are invalid for the distribution,
   no parameters are set.
   Notice that optional parameters are (re-)set to their default values if 
   not given for UNURAN standard distributions.

   @emph{Important:} Integer parameter must be given as doubles.
*/

/*
  =END
*/

/*---------------------------------------------------------------------------*/



