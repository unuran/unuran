/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_debug.h                                                         *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for debugging routines.    *
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

#include <stdio.h>

/*---------------------------------------------------------------------------*/
/* set debugging flag for generator                                          */
/*
 =INT
 */
int unur_set_debug( UNUR_PAR *parameters, unsigned debug );
int unur_chg_debug( UNUR_GEN *generator, unsigned debug );
int unur_set_default_debug( unsigned debug );
/*
 =END
 */

/* common debug flags                                                        */
#define UNUR_DEBUG_INIT    0x00000001u    /* bit  01 ... pameters of generator */
#define UNUR_DEBUG_SETUP   0x00000fffu    /* bits 02-12 ... setup            */
#define UNUR_DEBUG_ADAPT   0x00fff000u    /* bits 13-24 ... adaptive steps   */
#define UNUR_DEBUG_SAMPLE  0xff000000u    /* bits 25-32 ... trace sampling   */

#define UNUR_DEBUG_OFF     (0u)       /* switch off debugging information    */    
#define UNUR_DEBUG_ALL     (~0u)      /* write all avaivable information     */

/*---------------------------------------------------------------------------*/
/* global variable used to record errors                                     */
extern unsigned unur_errno;

/*---------------------------------------------------------------------------*/
/* manipulate output stream                                                  */
FILE *unur_set_stream( FILE *new_stream );
FILE *unur_get_stream( void );

/*---------------------------------------------------------------------------*/
/* warnings and error messages for given error number                        */
const char *unur_get_strerror ( const int unur_errno );

/*---------------------------------------------------------------------------*/




