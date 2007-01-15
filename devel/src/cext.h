/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cext.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method CEXT                               *
 *         (wrapper for Continuous EXTernal generators)                      *
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
   =expMETHOD  CEXT   wrapper for Continuous EXTernal generators

   =UP  Methods_for_CONT

   =REQUIRED 

   =SPEED depends on external generator

   =REINIT supported

   =DESCRIPTION
      CEXT is a wrapper for external generators for continuous
      univariate distributions. It allows the usage of external
      random variate generators within the UNURAN framework.

   =HOWTOUSE
      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_cext_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for new generator. 
*/

/*...........................................................................*/

int unur_cext_set_init( UNUR_PAR *parameters, int (*init)(UNUR_GEN *gen) );
/*
   Set initialization routine for external generator.
   
   @emph{Important:} The routine @var{init} must return
   @code{UNUR_SUCCESS} when the generator was initialized successfully
   and @code{UNUR_FAILURE} otherwise.
*/

int unur_cext_set_sample( UNUR_PAR *parameters, double (*sample)(UNUR_GEN *gen) );
/*
   Set sampling routine for external generator.
*/

void *unur_cext_get_params( UNUR_GEN *generator, size_t size );
/*
   Get pointer to memory block for storing parameters of external
   generator. A memory block of size @var{size} is automatically (re-)
   allocated if necessary and the pointer to this block is stored in
   the @var{generator} object. If one only needs the pointer to this
   memory block set @var{size} to @code{0}.

   Notice, that @var{size} is the size of the memory block and not the
   length of an array.

   @emph{Important:} This rountine should only be used in the
   initialization and sampling routine of the external generator. 
*/

double *unur_cext_get_distrparams( UNUR_GEN *generator );
/* */

int unur_cext_get_ndistrparams( UNUR_GEN *generator );
/*
   Get size of and pointer to array of parameters of underlying
   distribution in @var{generator} object.

   @emph{Important:} This rountine should only be used in the
   initialization and sampling routine of the external generator. 
*/

/*...........................................................................*/


/* =END */
/*---------------------------------------------------------------------------*/
