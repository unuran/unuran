/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_cstd.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method CSTD       *
 *         (generators for Continuous STanDard distributions)                *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only used in unur_methods.h                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#define UNUR_MAX_DIST_PARAMS  5    /* maximal numbers of parameters for distributions */

/*---------------------------------------------------------------------------*/
/* Information for constructing the generator                                */

struct unur_cstd_par { 
  char *definition; /* description of standard distribution */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_cstd_gen { 
  double dist_param[UNUR_MAX_DIST_PARAMS]; /* parameters for standard distribution */
  int    n_dist_param;  /* number of parameters for distribution             */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_cstd_new( char *definition );
/* get default parameters for generator                                      */

struct unur_gen *unur_cstd_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_cstd_sample( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_cstd_free( struct unur_gen *generator);
/* destroy generator object                                                  */


/** TODO: use unur_set_*** interface **/
/*  void cstd_set_1par(struct unur_par *par,double inpval); */
/*  void cstd_set_2par(struct unur_par *par,double inpval1,double inpval2); */

/*---------------------------------------------------------------------------*/

