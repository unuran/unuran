/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_dstd.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method DSTD       *
 *         (generators for Continuous STanDard distributions)                *
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

#define UNUR_MAX_DIST_PARAMS  5    /* maximal numbers of parameters for distributions */

/*---------------------------------------------------------------------------*/
/* Information for constructing the generator                                */

struct unur_dstd_par { 
  const char *sample_routine_name; /* name of sampling routine               */
  int  is_inversion;      /* indicate whether method is inversion method     */     
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_dstd_gen { 
  double *gen_param;      /* parameters for the generator                    */
  int     n_gen_param;    /* number of parameters for the generator          */

  int    *gen_iparam;     /* integer parameters for generator                */
  int     n_gen_iparam;   /* number of integer parameters for the generator  */

  double  umin;           /* cdf at left boundary of domain                  */
  double  umax;           /* cdf at right boundary of domain                 */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_dstd_new( struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_dstd_init( struct unur_par *parameters );
/* initialize new generator                                                  */

/** 
    double unur_dstd_sample( struct unur_gen *gen );
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

void unur_dstd_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_dstd_set_variant( struct unur_par *par, unsigned variant );
/* set variant of method                                                     */

int unur_dstd_chg_param( struct unur_gen *gen, double *params, int n_params );
/* change array of parameters for distribution                               */

#define unur_dstd_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

