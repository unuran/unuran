/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_stdgen.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines list of variants for special generators                   *
 *         for standard distributions.                                       *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in ../methods/cstd.c !!!                            *
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
#ifndef __UNUR_STDGEN_VARIANTS_H_SEEN
#define __UNUR_STDGEN_VARIANTS_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>
#include <unur_methods.h>

/*---------------------------------------------------------------------------*/

/* declare type of routine for sampling 
   from univariate continuous distribution                                   */
typedef double (_unur_sampling_routine_cont)( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/

/* the first entry [0] is the default variant */

/*---------------------------------------------------------------------------*/
/* exponential distribution                                                  */

#define CSTD_EXPONENTIAL_N_VAR  1

static _unur_sampling_routine_cont *
_cstd_exponential_var[CSTD_EXPONENTIAL_N_VAR] = {
  /* 0 */  unur_cstd_sample_exponential_inv              /* Inversion method */
};

#if UNUR_DEBUG & UNUR_DB_INFO
static char *_cstd_exponential_varname[CSTD_EXPONENTIAL_N_VAR] = {
  /* 0 */  "unur_cstd_sample_exponential_inv"            /* Inversion method */
};
#endif


/*---------------------------------------------------------------------------*/
/* gamma distribution                                                        */

#define CSTD_GAMMA_N_VAR  1

static _unur_sampling_routine_cont *
_cstd_gamma_var[CSTD_GAMMA_N_VAR] = {
  /* 0 */  unur_cstd_sample_gamma_gammarand              /* ???? method      */
};

#if UNUR_DEBUG & UNUR_DB_INFO
static char *_cstd_gamma_varname[CSTD_GAMMA_N_VAR] = {
  /* 0 */  "unur_cstd_sample_gamma_gammarand"            /* ???? method      */
};
#endif


/*---------------------------------------------------------------------------*/
/* normal distribution                                                       */

#define CSTD_NORMAL_N_VAR  1

static _unur_sampling_routine_cont *
_cstd_normal_var[CSTD_NORMAL_N_VAR] = {
  /* 0 */  unur_cstd_sample_normal_bm                   /* Box-Muller method */
};

#if UNUR_DEBUG & UNUR_DB_INFO
static char *_cstd_normal_varname[CSTD_NORMAL_N_VAR] = {
  /* 0 */  "unur_cstd_sample_normal_bm"                 /* Box-Muller method */
};
#endif


/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_STDGEN_VARIANTS_H_SEEN */
/*---------------------------------------------------------------------------*/
