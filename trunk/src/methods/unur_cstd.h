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

struct unur_cstd_par { 
  unsigned int variant;   /* indicates generator to be used                  */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_cstd_gen { 
  double *pdf_param;      /* parameters for standard distribution            */
  int     n_pdf_param;    /* number of parameters for distribution           */
  unsigned int variant;   /* indicates generator to be used                  */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_cstd_new( struct unur_distr *distr );
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

