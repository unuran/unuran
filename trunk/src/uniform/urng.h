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
#ifndef URNG_H_SEEN
#define URNG_H_SEEN
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
      @ @ @ double (*getrand)(void *state);
      @ @ @ void *state;
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
        urng->state = my_parameters;                                   
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

/*---------------------------------------------------------------------------*/


/* =END */








/*---------------------------------------------------------------------------*/
/* Unified interface for uniform random number generators.                   */
/* EXPERIMENTAL!                                                             */


double unur_urng_sample (UNUR_URNG *urng);
/*
   Get a uniform random number from @var{urng}.
   If the NULL pointer is given, the default uniform generator is
   used. 
*/


int unur_urng_reset (UNUR_URNG *urng);
/*
   Reset @var{urng} object. 
   The routine tries two ways to reset the generator (in this order):

   (1) It uses the reset function given by an unur_urng_set_reset()
       call. 

   (2) It uses the seed given by the last unur_urng_seed() call (which
       requires a seeding function given by a unur_urng_set_seed()
       call). If no seed was given explicitly (or implicitly for some
       kinds of uniform random number generators) the generator cannot
       be reset.

   If neither of the two methods work resetting of the generator is
   not possible and an error code is returned.

   If the NULL pointer is given, the default uniform generator is
   reset.
*/



/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

/* generic API to create a new URNG (uniform random number generator) object */

UNUR_URNG *unur_urng_new( double (*sampleunif)(void *state), void *state );
/*
   Get a new URNG object. 
   @var{sampleunif} is a function to the uniform sampling routine,
   @var{state} a pointer to its arguments which usually contains the
   state variables of the generator.

   Functions @var{sampleunif} with a different type for @var{p} or
   without an argument at all also work. A typecast might be necessary
   to avoid compiler warnings.

   Functions @var{sampleunif} does not have any argument use NULL for
   @var{state}.
   
   @emph{Important:} @var{sampleunif} must not be the NULL pointer.

   There are appropriate calls that simplifies the task of creating
   URNG objects for some libraries with uniform random number
   generators.
*/

int unur_urng_set_seed( UNUR_URNG *urng, void (*setseed)(void *state, unsigned long seed) );
/*
   Set function to seed generator @var{urng} (if available).
*/

int unur_urng_set_anti( UNUR_URNG *urng, void (*setanti)(void *state, int anti) );
/*
   Set function to switch the antithetic flag of generator @var{urng}
   (if available).
*/

int unur_urng_set_reset( UNUR_URNG *urng, void (*reset)(void *state) );
/* 
   Set function for reseting the uniform random number generator
   @var{urng} (if available).
*/

int unur_urng_set_nextsub( UNUR_URNG *urng, void (*nextsub)(void *state) );
/*
   Set function that allows jumping to start of the next substream of
   @var{urng} (if available).
*/

int unur_urng_set_resetsub( UNUR_URNG *urng, void (*resetsub)(void *state) );
/*
   Set function that allows jumping to start of the current substream
   of @var{urng} (if available).
*/

int unur_urng_set_delete( UNUR_URNG *urng, void (*delete)(void *state) );
/*
   Set function for destroying @var{urng} (if available).
*/

/*---------------------------------------------------------------------------*/

/* generic API to handle URNG (uniform random number generator) object */

int unur_urng_seed (UNUR_URNG *urng, unsigned long seed);
/*
   Set @var{seed} for generator @var{urng} (reseed)
   (if available).
   It returns an error code if this is not possible. 

   If the NULL pointer is given, the default uniform generator is
   seeded (if possible).

   @emph{Notice}: Seeding should be done only once for a particular
   generator (except for resetting it to the initial state).
   Expertise is required when seeding is used to get independent
   streams. Thus we recommend appropriate libraries for this task,
   e.g. Pierre L'Ecuyer's @file{RngStreams} package. For this library
   only a package seed can be set and thus the unur_urng_seed() call
   will not have any effect to generators of this type. Use
   unur_urng_reset() or unur_urng_rngstream_new() instead, depending
   whether one wants to reset the stream or get a new stream that is
   independent from the previous ones.
*/

int unur_urng_anti (UNUR_URNG *urng, int anti);
/*
   Switch to antithetic random numbers in @var{urng}.
   (if available).
   It returns an error code if this is not possible. 

   If the NULL pointer is given, the default uniform generator is
   reset (if possible).
*/

int unur_urng_nextsub (UNUR_URNG *urng);
/*
   Jump to start of next substream of @var{urng}.
   It returns an error code if this is not possible. 

   If the NULL pointer is given, the default uniform generator is set
   to the start of the next substream (if possible).
*/

int unur_urng_resetsub (UNUR_URNG *urng);
/*
   Jump to start of the current substream of @var{urng}.
   It returns an error code if this is not possible. 

   If the NULL pointer is given, the default uniform generator is set
   to the start of the current substream (if possible).
*/

int unur_urng_free (UNUR_URNG *urng);
/* 
   Destroy @var{urng} object.
   It returns an error code if this is not possible. 

   If the NULL is given, this function does nothing as it is not
   a good idea to destroy the default uniform generator.

   @emph{Warning:} This call must be used with care. The @var{urng}
   object must not be used by any existing generator object!
*/

int unur_gen_seed (UNUR_GEN *gen, unsigned long seed);
/* */

int unur_gen_anti (UNUR_GEN *gen, int anti);
/* */

int unur_gen_reset (UNUR_GEN *gen);
/* */

int unur_gen_nextsub (UNUR_GEN *gen);
/* */

int unur_gen_resetsub (UNUR_GEN *gen);
/* 
   Analogous to unur_urng_seed(), unur_urng_anti(), unur_urng_reset(),
   unur_nextsub(), and unur_urng_resetsub(). 
   However, these act on the URNG object used by generator object
   @var{gen}. 

   @emph{Warning:} These three calls should be used with care as it
   influences all generator objects that share the same URNG object!
*/

/*---------------------------------------------------------------------------*/

/* Interface to particular uniform RNG libraries */


UNUR_URNG *unur_urng_fvoid_new( double (*random)(void), int (*reset)(void) );
/*
   Make a URNG object for a genertor that consists of a single
   function call with a global state variable.

   @emph{Notice:} If independent versions of the same URNG should be
   used, copies of the subroutine with different names has to be
   implement in the program code.

   If there is no reset function use NULL for the second argument.
   
   UNURAN contains some build-in URNGs of this type in directory
   @file{src/uniform/}.
*/

#ifdef UNURAN_HAS_PRNG

UNUR_URNG *unur_urng_prng_new( const char *prngstr );
/*
   Make object for URNGs from Otmar Lendl's @file{prng} package. 
   @var{prngstr} is a string that contains the necessary information
   to create a uniform random number generator. For the format of this
   string see the @file{prng} user manual.

   The @file{prng} library provides a very flexible way to sample form
   arbitrary URNGs by means of an object oriented programing
   paradigma. Similarly to the UNURAN library independent generator
   objects can be build and used. The library has been developed by
   and implemented by Otmar Lendl as member of the pLab group at the
   university of Salzburg (Austria, EU). 

   It is available via anonymous ftp from
   @uref{http://statistik.wu-wien.ac.at/prng/}
   or from the pLab site at
   @uref{http://random.mat.sbg.ac.at/}.
*/

UNUR_URNG *unur_urng_prngptr_new( struct prng *urng );
/*
   Similar to unur_urng_prng_new() except it uses a pointer to a
   generator object as returned by @code{prng_new(prngstr)};
   see @file{prng} manual for details.
*/

#endif



#ifdef UNURAN_HAS_RNGSTREAMS

UNUR_URNG *unur_urng_rngstream_new( const char *urngstr );
/*
   Make object for URNGs from Pierre L'Ecuyer's @file{RngStream}
   library. This library provides multiple independent streams of
   pseudo-random numbers and is available from
   @uref{http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/}.
   @var{urngstr} is an arbitrary string to label a stream. It need not
   be unique.
*/

UNUR_URNG *unur_urng_rngstreamptr_new( RngStream rngstream );
/*
   Similar to unur_urng_rngstream_new() except it uses a pointer to a 
   generator object as returned by @code{RngStream_CreateStream()}.
*/

#endif



#ifdef UNURAN_HAS_GSL

UNUR_URNG *unur_urng_gsl_new( const gsl_rng_type *urngtype );
/*
   Make object for URNGs from the @file{GSL} (GNU Scientific Library).
   @var{urng} is the type of the chosen generator as described in the
   GSL manual. This library is available from
   @uref{http://www.gnu.org/software/gsl/}.
*/

UNUR_URNG *unur_urng_gslptr_new( gsl_rng *urng );
/*
   Similar to unur_urng_gsl_new() except it uses a pointer to a
   generator object as returned by @code{gsl_rng_alloc(rng_type)};
   see @file{GSL} manual for details.

   @emph{Notice}: There is a subtle but important difference between
   these two calls. When a generator object is created by a 
   unur_urng_gsl_new() call, then resetting of the generator works.
   When a generator object is created by a unur_urng_gslptr_new()
   call, then resetting only works after a
   @code{unur_urng_seed(urng,myseed)} call. 
*/

#endif

/*---------------------------------------------------------------------------*/
#endif   /* #if UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#endif  /* URNG_H_SEEN */
/*---------------------------------------------------------------------------*/


