/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_misc.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines function prototypes for miscelleanous routines            *
 *         parameters in generator objects.                                  *
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
#ifndef __UNUR_MISC_H_SEEN
#define __UNUR_MISC_H_SEEN
/*---------------------------------------------------------------------------*/

#include <stdio.h>

/*---------------------------------------------------------------------------*/
/* Defining infinity                                                         */
/* (we use the largest possible value to indicate infinity)                  */
#include <math.h>
#define UNUR_INFINITY  HUGE_VAL     

/*---------------------------------------------------------------------------*/
/* set debugging flag for generator                                          */
int unur_set_debug( struct unur_par *parameter, unsigned debug );
int unur_set_default_debug( unsigned debug );
extern unsigned _unur_default_debugflag;     /* default debugging flags      */

/* common debug flags                                                        */
#define UNUR_DEBUG_INIT    0x00000001u    /* bit  01 ... pameters of generator */
#define UNUR_DEBUG_SETUP   0x00000fffu    /* bits 02-12 ... setup            */
#define UNUR_DEBUG_ADAPT   0x00fff000u    /* bits 13-24 ... adaptive steps   */
#define UNUR_DEBUG_SAMPLE  0xff000000u    /* bits 25-32 ... trace sampling   */

#define UNUR_DEBUG_OFF     (0u)       /* switch off debugging information    */    
#define UNUR_DEBUG_ALL     (~0u)      /* write all avaivable information     */

/*---------------------------------------------------------------------------*/
/* manipulate output stream                                                  */
FILE *unur_set_stream( FILE *new_stream );
FILE *unur_get_stream( void );

/*---------------------------------------------------------------------------*/
/* warnings and error messages for given error number                        */
const char *unur_get_strerror ( const int unur_errno );

/*---------------------------------------------------------------------------*/
/* set, get or change uniform RNG for generator                              */

int unur_set_urng( struct unur_par *par, UNUR_URNG_TYPE urng );
UNUR_URNG_TYPE unur_chg_urng( struct unur_gen *gen, UNUR_URNG_TYPE urng );
UNUR_URNG_TYPE unur_get_urng( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* get and set default uniform RNG                                           */
/* (defined in src/utils/urng.c)                                             */

UNUR_URNG_TYPE unur_get_default_urng( void );
UNUR_URNG_TYPE unur_set_default_urng( UNUR_URNG_TYPE urng_new );

/*---------------------------------------------------------------------------*/
/* get dimension of generator for (multivariate) distribution                */

int unur_get_dimension( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* get type of transformation method                                         */

#define unur_is_discr(gen) ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ? 1 : 0 )
#define unur_is_cont(gen)  ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_CONT)  ? 1 : 0 )
#define unur_is_vec(gen)   ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_VEC)   ? 1 : 0 )

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_MISC_H_SEEN */
/*---------------------------------------------------------------------------*/








