/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_stdgen.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for special generators     *
 *         for standard distributions.                                       *
 *                                                                           *
 *   USAGE:                                                                  *
 *         included in all generator source files.                           *
 *         included in ../methods/cstd.c.                                    *
 *         required for every application of special generators when not     *
 *         invoked via the normal UNURAN interface.                          *
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
#ifndef __UNUR_STDGEN_H_SEEN
#define __UNUR_STDGEN_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>
#include <unur_urng.h>

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* normal distribution                                                       */
double nbm(UNUR_URNG_TYPE urng);

/*---------------------------------------------------------------------------*/
/* gamma distribution                                                        */
double gammarand(double a, UNUR_URNG_TYPE urng);


/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_STDGEN_H_SEEN */
/*---------------------------------------------------------------------------*/
