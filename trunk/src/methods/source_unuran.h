/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_unuran.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and declares structures and function prototypes    *
 *         for all UNURAN source files                                       *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source files.                                    *
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
#ifndef UNURAN_SOURCE_H_SEEN
#define UNURAN_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* include main header files                                                 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#  include <config.h>
#else
#  error "config.h" required
#endif

#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif

#include <in_unuran.h>

#include <parser_source.h>
#include <functparser_source.h>

#include <source_struct.h>
#include <source_cookies.h>
#include <source_fp.h>
#include <source_methods.h>

#include <distr_source.h>
#include <x_debug_source.h>
#include <x_gen_source.h>
#include <x_math_source.h>
#include <x_misc_source.h>
#include <x_umalloc_source.h>

#ifdef WITH_DMALLOC
#  include <dmalloc.h>
#endif

/*---------------------------------------------------------------------------*/
#endif  /* end UNURAN_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
