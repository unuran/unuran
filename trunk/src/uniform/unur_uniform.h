/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_uniform.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for included uniform random number generators *
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
#ifndef UNUR_UNIFORM_H_SEEN
#define UNUR_UNIFORM_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unuran_config.h>

/*---------------------------------------------------------------------------*/

/* Combined multiple recursive generator by Pierre L'Ecuyer and Renee Touzin */
/* Copyright (c) 2002 Renee Touzin.                                          */

double unur_urng_MRG31k3p (void);
int unur_urng_MRG31k3p_seed (long seed);
int unur_urng_MRG31k3p_reset (void);

/* Linear congruential generator by Fishman and Moore                        */
/* m = 2^31-1, a = 742938285, c = 0.                                         */

double unur_urng_fish (void);
int unur_urng_fish_seed (long seed);
int unur_urng_fish_reset (void);

/* Linear congruential generator "Minimal Standard"                          */
/* m = 2^31-1, a = 16807, c = 0.                                             */

double unur_urng_mstd (void);
int unur_urng_mstd_seed (long seed);
int unur_urng_mstd_reset (void);

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_UNIFORM_H_SEEN */
/*---------------------------------------------------------------------------*/
