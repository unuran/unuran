/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_distribution.h                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros, typedefs and constants and for PDF, CDF, etc.     *
 *         of distribtions.                                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in distribution source files.                       *
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
#ifndef __SOURCE_DISTRIBUTIONS_LIB_H_SEEN
#define __SOURCE_DISTRIBUTIONS_LIB_H_SEEN
/*---------------------------------------------------------------------------*/

#include <source_unuran.h>
#include <source_specfunct.h>

/*---------------------------------------------------------------------------*/
/* Macros                                                                    */

/*---------------------------------------------------------------------------*/
/** TODO **/
#define NOT_UNIMODAL  0
#define RETURN_NULL   0
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* set routine for sampling                                                  */
#define _unur_cstd_set_sampling_routine(par,gen,routine) \
   do { \
     if ((gen)==NULL) return 1;                    /* test existence only  */ \
     (gen)->sample.cont = (routine);                 /* set pointer        */ \
     if (par) (par)->data.cstd.sample_routine_name = #routine;  /* set routine name */ \
   } while (0)

#define _unur_dstd_set_sampling_routine(par,gen,routine) \
   do { \
     if ((gen)==NULL) return 1;                    /* test existence only  */ \
     (gen)->sample.discr = (routine);                /* set pointer        */ \
     if (par) (par)->data.dstd.sample_routine_name = #routine;  /* set routine name */ \
   } while (0)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_DISTRIBUTIONS_LIB_H_SEEN */
/*---------------------------------------------------------------------------*/
