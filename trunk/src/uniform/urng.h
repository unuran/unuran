/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros and function prototypes for using uniform         *
 *         random number generators inside UNURAN.                           *
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

/*---------------------------------------------------------------------------*/
#ifndef X_URNG_H_SEEN
#define X_URNG_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unuran_config.h>

/* 
   =NODE  URNG  Using uniform random number generators

   =UP TOP [50]

   =DESCRIPTION
      Each generator has a pointer to a uniform (pseudo-) random number
      generator (URNG). It can be set via the unur_set_urng() call. It is
      also possible to read this pointer via unur_get_urng() or change the
      URNG for an existing generator object by means of unur_chg_urng();
      
      By this very flexible concept it is possible that each generator has
      its own (independent) URNG or several generators can share the same
      URNG. 
      
      If no URNG is provided for a parameter or generator object a default
      generator is used which is the same for all generators. This URNG is
      defined in @file{unuran_config.h} at compile time. A pointer to
      this default URNG can be obtained via unur_get_default_urng().
      Nevertheless, it is also possible to overwrite this default URNG by
      another one by means of the unur_set_default_urng() call. However,
      this only takes effect for new parameter objects.
      
      The pointer to a URNG is of type @code{UNUR_URNG*}. Its definition 
      depends on the compilation switch @code{UNUR_URNG_TYPE} in @
      @file{unuran_config.h}. 
      Currently we have two possible switches (other values would result
      in a compilation error):
      
      @enumerate

      @item
      UNUR_URNG_TYPE == UNUR_URNG_FVOID
      
      This uses URNGs of type @code{double uniform(void)}.
      If independent versions of the same URNG should be used, a copy of
      the subroutine has to be implement in the program code (with
      different names, of course).
      UNURAN contains some build-in URNGs of this type in directory
      @file{src/uniform/}.
      
      @item
      UNUR_URNG_TYPE == UNUR_URNG_PRNG
      
      This uses the URNGs from the @code{prng} library. It provides a very
      flexible way to sample form arbitrary URNGs by means of an object
      oriented programing paradigma. Similarly to the UNURAN library
      independent generator objects can be build and used.
      Here @code{UNUR_URNG*} is simply a pointer to such a uniform
      generator object.
      
      This library has been developed by the pLab group at the university
      of Salzburg (Austria, EU) and implemented by Otmar Lendl.
      It is available via anonymous ftp from
      @uref{http://statistik.wu-wien.ac.at/prng/}
      or from the pLab site at
      @uref{http://random.mat.sbg.ac.at/}.
      
      @item
      UNUR_URNG_TYPE == UNUR_URNG_RNGSTREAM

      Use Pierre L'Ecuyer's @code{RngStream} library for multiple 
      independent streams of pseudo-random numbers. 
      It is available from
      @uref{http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/}.        

      @item
      UNUR_URNG_TYPE == UNUR_URNG_GSL

      Use the URNG from the GNU Scientific Library (GSL).
      It is available from
      @uref{http://www.gnu.org/software/gsl/}.                               

      @item
      UNUR_URNG_TYPE == UNUR_URNG_GENERIC

      This a generic interface with limited support.
      It uses a structure to store both a function call of type
      @code{double urng(void*)} and a void pointer to the parameter list.
      Both pointers must be set directly using the structure
      @code{struct unur_urng_generic}
      (there are currently no calls that support this URNG type).
      It is defined as

      @smallexample
      struct unur_urng_generic @{
      @ @ @ double (*getrand)(void *params);
      @ @ @ void *params;
      @};
      @end smallexample

      All functions and parameters should be set at run time:
        @enumerate

	@item
	Allocate variable of type @code{struct unur_urng_generic}
	(or of type @code{UNUR_URNG}, which is the same):
	@smallexample
	UNUR_URNG *urng;                                                
	urng = malloc(sizeof(UNUR_URNG));
	@end smallexample

	@item
	Set function of type @code{double (*rand)(void *)} for
	sampling from URNG:
	@smallexample
        urng->getrand = my_uniform_rng;                                 
	@end smallexample

	@item
	Set pointer to parameters of for this function                    
	(or NULL if no parameters are required):                          
	@smallexample
        urng->params = my_parameters;                                   
	@end smallexample

	@item
	Use this URNG:
	@smallexample
        unur_urng_set_default(urng);             (set default generator)
        unur_urng_set_default_aux(urng);    (set default aux. generator)
	@end smallexample

	Notice that this must be done before UNURAN generator object      
	are created. Of course urng can also used just for a particular   
	generator object. Use the following and similar calls:
	@smallexample
        unur_set_urng(par,urng);                                        
        unur_chg_urng(gen,urng);                                        
	@end smallexample

        @end enumerate

      @end enumerate

      It is possible to use other interfaces to URNGs without much
      troubles. All changes have to be done in file 
      @file{unuran_config.h}. If you need such a new interface 
      please feel free to contact the authors of the UNURAN library.
      
      Some generating methods provide the possibility of correlation
      induction. To use this feature a second auxiliary URNG is required.
      It can be set and changed by the unur_set_urng_aux() and
      unur_chg_urng_aux() call, respectively. Since the auxiliary
      generator is by default the same as the main generator, the
      auxiliary URNG must be set after any unur_set_urng() or
      unur_chg_urng() call! Since in special cases mixing of two URNG
      might cause problems, we supply a default auxiliary generator that
      can be used by the unur_use_urng_aux_default() call (after the main
      URNG has been set). This default auxiliary generator can be changed
      with analogous calls as the (main) default uniform generator.

   =END

*/

/*---------------------------------------------------------------------------*/
/* set, get or change uniform RNG for generator                              */

/* =ROUTINES */

/*---------------------------------------------------------------------------*/
/* get and set default uniform RNG                                           */

/* ==DOC
   @subheading Default uniform RNGs
*/

UNUR_URNG *unur_get_default_urng( void );
/*
  Get the pointer to the default URNG. The default URNG is used by all
  generators where no URNG was set explicitly by a unur_set_urng()
  call.
*/

UNUR_URNG *unur_set_default_urng( UNUR_URNG *urng_new );
/*
  Change the default URNG for new parameter objects. 
*/


UNUR_URNG *unur_set_default_urng_aux( UNUR_URNG *urng_new );
/* */

UNUR_URNG *unur_get_default_urng_aux( void );
/*
  Analogous calls for default auxiliary generator.
*/

/*---------------------------------------------------------------------------*/

/* ==DOC
   @subheading Uniform RNGs for generator objects
*/

int unur_set_urng( UNUR_PAR *parameters, UNUR_URNG *urng );
/*
  Use the URNG @code{urng} for the new generator. This overwrite the
  default URNG. It also sets the auxiliary URNG to @code{urng}.

  @emph{Important}: For multivariate distributions that use 
  marginal distributions this call does not work properly.
  It is then better first to create the generator object (by
  a unur_init() call) and then change the URNG by means of 
  unur_chg_urng().
*/

UNUR_URNG *unur_chg_urng( UNUR_GEN *generator, UNUR_URNG *urng );
/*
  Change the URNG for the given generator. It returns the pointer to
  the old URNG that has been used by the generator.
  It also changes the auxiliary URNG to @code{urng} and thus
  overwrite the last unur_chg_urng_aux() call.
*/

UNUR_URNG *unur_get_urng( UNUR_GEN *generator );
/*
  Get the pointer to the URNG that is used by the generator.
  This is usefull if two generators should share the same URNG.
*/

int unur_set_urng_aux( UNUR_PAR *parameters, UNUR_URNG *urng_aux );
/*
  Use the auxiliary URNG @code{urng_aux} for the new generator. 
  (Default is the default URNG or the URNG from the last
  unur_set_urng() call. Thus if the auxiliary generator should be
  different to the main URNG, unur_set_urng_aux() must be called after
  unur_set_urng(). 
  The auxiliary URNG is used as second stream of uniform random
  number for correlation induction.
  It is not possible to set an auxiliary URNG for a method that does
  not use one (i.e. the call returns an error code).
*/

int unur_use_urng_aux_default( UNUR_PAR *parameters );
/* 
   Use the default auxiliary URNG.
   (It must be set after unur_get_urng().)
   It is not possible to set an auxiliary URNG for a method that does
   not use one (i.e. the call returns an error code).
*/

int unur_chgto_urng_aux_default( UNUR_GEN *generator );
/*
   Switch to default auxiliary URNG.
   (It must be set after unur_get_urng().)
   It is not possible to set an auxiliary URNG for a method that does
   not use one (i.e. the call returns an error code).
*/

UNUR_URNG *unur_chg_urng_aux( UNUR_GEN *generator, UNUR_URNG *urng_aux );
/*
  Change the auxiliary URNG for the given generator. It returns the
  pointer to the old auxiliary URNG that has been used by the
  generator. It has to be called after each unur_chg_urng() when the 
  auxiliary URNG should be different from the main URNG.
  It is not possible to change the auxiliary URNG for a method that
  does not use one (i.e. the call NULL).
*/

UNUR_URNG *unur_get_urng_aux( UNUR_GEN *generator );
/*
  Get the pointer to the auxiliary URNG that is used by the
  generator. This is usefull if two generators should share the same
  URNG.
*/

/* =END */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* X_URNG_H_SEEN */
/*---------------------------------------------------------------------------*/
