/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_srou.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method SROU       *
 *         (Simple universal generator, ratio-of-uniforms method)            *
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
/* Information for constructing the generator                                */

struct unur_srou_par { 
  double Fmode;              /* cdf at mode                                  */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_srou_gen { 
  _UNUR_FUNCTION_CONT *pdf;  /* pointer to p.d.f.                            */
  double  mode;              /* location of mode                             */
  double *pdf_param;         /* parameters of the pdf                        */
  int     n_pdf_param;       /* number of parameters of the pdf              */

  double  um;                /* height of rectangle: square root of f(mode)  */
  double  vl, vr;            /* left and right boundary of rectangle         */
  double  xl, xr;            /* ratios vl/um and vr/um                       */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_srou_new( struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_srou_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_srou_sample( struct unur_gen *generator );
double unur_srou_sample_check( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_srou_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_srou_set_Fmode( struct unur_par *par, double Fmode );
/* set cdf at mode                                                           */

int unur_srou_set_verify( struct unur_par *par, int verify );
/* turn verifying of algorithm while sampling on/off                         */

int unur_srou_set_usesqueeze( struct unur_par *par, int usesqueeze );
/* set flag for using universal squeeze (default: off)                       */

int unur_srou_set_usemirror( struct unur_par *par, int usemirror );
/* set flag for using mirror principle (default: off)                        */

#define unur_srou_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

