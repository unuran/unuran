/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_distribution.h                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros, typedefs and constants and for p.d.f., c.d.f.,    *
 *         etc. of distribtions.                                             *
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

#include <unuran_config.h>
#include <source_stddistr.h>
#include <source_specfunct.h>
#include <unuran_distributions.h>

/*---------------------------------------------------------------------------*/
/* Macros                                                                    */

/*---------------------------------------------------------------------------*/
/** TODO **/
#define NOT_UNIMODAL  0
#define RETURN_NULL   0
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* (log of) normalization constant for p.d.f */
#define LOGNORMCONSTANT  params[UNUR_DISTR_MAXPARAMS]
#define NORMCONSTANT     params[UNUR_DISTR_MAXPARAMS]

/*---------------------------------------------------------------------------*/
/* set routine for sampling                                                  */
#define _unur_cstd_set_sampling_routine(par,gen,routine) \
   do { \
     if ((gen)==NULL) return 1;                    /* test existence only  */ \
     (gen)->sample.cont = (routine);                 /* set pointer        */ \
     (par)->data.cstd.sample_routine_name = #routine;  /* set routine name */ \
   } while (0)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_DISTRIBUTIONS_LIB_H_SEEN */
/*---------------------------------------------------------------------------*/
