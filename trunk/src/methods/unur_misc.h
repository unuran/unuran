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
 *         only used in unur_methods.h                                       *
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


/*---------------------------------------------------------------------------*/
/* set debugging flag for generator                                          */
int unur_set_debug( struct unur_par *parameter, unsigned debug );

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
#endif  /* __UNUR_MISC_H_SEEN */
/*---------------------------------------------------------------------------*/








