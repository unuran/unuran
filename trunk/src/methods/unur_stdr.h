/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_stdr.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method STDR       *
 *         (transformed density rejection with universal bounds)             *
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

struct unur_stdr_par { 
  double Fmode;              /* cdf at mode                                  */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_stdr_gen { 
  _UNUR_FUNCTION_CONT *pdf;  /* pointer to p.d.f.                            */
  double  mode;              /* location of mode                             */
  double *pdf_param;         /* parameters of the pdf                        */
  int     n_pdf_param;       /* number of parameters of the pdf              */

  double  fm;                /* pdf at mode                                  */
  double  um;                /* sqrt of pdf at mode                          */
  double  vl, vr;            /* parameters for hat function                  */
  double  xl, xr;            /* partition points of hat                      */
  double  al, ar;            /* areas below hat in first and secont part     */
  double  A;                 /* area below hat                               */
  double  Aleft, Ain;        /* areas below hat in left tails and inside domain of pdf */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_stdr_new( struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_stdr_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_stdr_sample( struct unur_gen *generator );
double unur_stdr_sample_check( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_stdr_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_stdr_set_Fmode( struct unur_par *par, double Fmode );
/* set cdf at mode                                                           */

int unur_stdr_set_verify( struct unur_par *par, int verify );
/* turn verifying of algorithm while sampling on/off                         */

int unur_stdr_set_usesqueeze( struct unur_par *par, int usesqueeze );
/* set flag for using universal squeeze (default: off)                       */

#define unur_stdr_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

