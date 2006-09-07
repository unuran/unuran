/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: debug.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for debugging routines.    *
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

/* 
   =NODE  Debug  Debugging

   =UP TOP [70]

   =DESCRIPTION
      The UNURAN library has several debugging levels which can be
      switched on/off by debugging flags. This debugging feature can be
      enabled by defining the macro @code{UNUR_ENABLE_LOGGING} in 
      @file{unuran_config.h}.
      The debugging levels range from print a short description of the build
      generator object to a detailed description of hat functions and
      tracing the sampling routines. The output is print onto the output
      stream obtained by unur_get_stream() (see also ?).
      These flags can be set or changed by the respective calls
      unur_set_debug() and unur_chg_debug() independently for each generator. 
      The default debugging flags are given by the macro
      @code{UNUR_DEBUGFLAG_DEFAULT} in @file{unuran_config.h}.
      This default can be overwritten at run time by a
      unur_set_default_debug() call.
      
      Off course these debugging flags
      depend on the chosen method. Since most of these are merely for
      debugging the library itself, a description of the flags are given
      in the corresponding source files of the method.
      Nevertheless, the following flags can be used with all methods.


      Common debug flags:

      @ftable @code
      @item UNUR_DEBUG_OFF
      @findex UNUR_DEBUG_OFF
      switch off all debuging information
      @item UNUR_DEBUG_ALL
      @findex UNUR_DEBUG_ALL
      all avaivable information
      @item UNUR_DEBUG_INIT
      @findex UNUR_DEBUG_INIT
      parameters of generator object after initialization
      @item UNUR_DEBUG_SETUP
      @findex UNUR_DEBUG_SETUP
      data created at setup
      @item UNUR_DEBUG_ADAPT
      @findex UNUR_DEBUG_ADAPT
      data created during adaptive steps 
      @item UNUR_DEBUG_SAMPLE
      @findex UNUR_DEBUG_SAMPLE
      trace sampling   
      @end ftable


      Almost all routines check a given pointer before they read from or write
      to the given address. This does not hold for time-critical routines
      like all sampling routines. Then you are responsible for checking a
      pointer that is returned from a unur_init() call.
      However it is possible to turn on checking for invalid NULL pointers
      even in such time-critical routines by defining
      @code{UNUR_ENABLE_CHECKNULL} in @file{unuran_config.h}.

      Another debugging tool used in the library are magic cookies that
      validate a given pointer. It produces an error whenever a given
      pointer points to an object that is invalid in the context.
      The usage of magic cookies can be switched on by defining
      @code{UNUR_COOKIES} in @file{unuran_config.h}.

   =END
*/

/* common debug flags                                                        */
#define UNUR_DEBUG_OFF     (0u)       /* switch off debugging information    */    
#define UNUR_DEBUG_ALL     (~0u)      /* write all avaivable information     */

#define UNUR_DEBUG_INIT    0x00000001u    /* bit  01 ... parameters of gen.  */
#define UNUR_DEBUG_SETUP   0x00000fffu    /* bits 02-12 ... setup            */
#define UNUR_DEBUG_ADAPT   0x00fff000u    /* bits 13-24 ... adaptive steps   */
#define UNUR_DEBUG_SAMPLE  0xff000000u    /* bits 25-32 ... trace sampling   */

/*---------------------------------------------------------------------------*/
/* set debugging flag for generator                                          */

/* =ROUTINES */

int unur_set_debug( UNUR_PAR *parameters, unsigned debug );
/*
  Set debugging flags for generator.
*/

int unur_chg_debug( UNUR_GEN *generator, unsigned debug );
/*
  Change debugging flags for generator.
*/

int unur_set_default_debug( unsigned debug );
/*
  Overwrite the default debugging flag.
*/

/* =END */

/*---------------------------------------------------------------------------*/




