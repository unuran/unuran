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
#include <unur_methods.h>

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* exponential distribution                                                  */

double unur_cstd_sample_exponential_inv( struct unur_gen *gen );
/* Inversion method                                                          */

/*---------------------------------------------------------------------------*/
/* gamma distribution                                                        */
double unur_cstd_sample_gamma_gammarand( struct unur_gen *gen );
/* ???? method                                                    */

/*---------------------------------------------------------------------------*/
/* normal distribution                                                       */

double unur_cstd_sample_normal_bm( struct unur_gen *gen );
/* Box-Muller method                                                         */

double unur_cstd_sample_normal_pol( struct unur_gen *gen );
/* Polarmethod with rejection                                                */

double unur_cstd_sample_normal_quo( struct unur_gen *gen );
/* Ratio-of-uniforms method with squeeze                                     */

double unur_cstd_sample_normal_nquo( struct unur_gen *gen );
/* "Naive" ratio-of-uniforms method                                          */

double unur_cstd_sample_normal_leva( struct unur_gen *gen );
/* Ratio-of-uniforms method  with quadratic bounding curves                  */

double unur_cstd_sample_normal_kr( struct unur_gen *gen );
/* Kindermann-Ramage method                                                  */

double unur_cstd_sample_normal_acr( struct unur_gen *gen );
/* Acceptance-complement ratio                                               */

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_STDGEN_H_SEEN */
/*---------------------------------------------------------------------------*/
