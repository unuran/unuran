/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      unur_set.h                                                   *
 *                                                                           *
 *****************************************************************************
 * $Id$                                                                      *
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
#ifndef __UNUR_DEFS_H_SEEN
#define __UNUR_DEFS_H_SEEN
/*---------------------------------------------------------------------------*/

#include <math.h>

/*---------------------------------------------------------------------------*/
/* Defining infinity                                                         */
/* (we use the largest possible value to indicate infinity)                  */

#define UNUR_INFINITY  HUGE_VAL     

/*---------------------------------------------------------------------------*/
/* invoke uniform random number generator                                    */

#define UNUR_URNG_LINKED   1     /* link routine at compilation              */
#define UNUR_URNG_POINTER  2     /* use a pointer to a routine               */
#define UNUR_URNG_PRNG     3     /* use a pointer to gen. from prng-2.2      */

/*---------------------------------------------------------------------------*/
/* debugging flags                                                           */

#define UNUR_DB_STDERR     0x001UL  /* write warnings and errors on stderr   */
#define UNUR_DB_LOG        0x002UL  /* write warnings and infos into logfile */
#define UNUR_DB_COOKIES    0x004UL  /* use magic cookies                     */
#define UNUR_DB_CHECKNULL  0x008UL  /* check for NULL pointer                */
#define UNUR_DB_CHECKARGS  0x010UL  /* check arguments                       */
#define UNUR_DB_INFO       0x020UL  /* write info about generator into logfile */

/*---------------------------------------------------------------------------*/

#include <unuran_config.h>

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_DEFS_H_SEEN */
/*---------------------------------------------------------------------------*/

