/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_utdr.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method UTDR       *
 *         (Universal Transformed Density Rejection; 3-point method)         *
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
/* Information for constructing the generator                                */

struct unur_utdr_par { 

  double (*pdf)(double x, double *pdf_param, int n_pdf_param);  /* pointer to pdf */
  double  pdf_area;     /* approximate area below the pdf                    */
  double  mode;         /* location of mode                                  */
  double  il;           /* left border of the domain                         */
  double  ir;           /* right border of the domain                        */
  double *pdf_param;    /* parameters of the pdf                             */
  int     n_pdf_param;  /* number of parameters of the pdf                   */
  double  c_factor;     /* constant for choosing the design points           */
  double  delta_factor; /* constant for choosing delta to replace the tangent*/
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_utdr_gen { 

  double (*pdf)(double x, double *pdf_param, int n_pdf_param);  /* pointer to pdf */
  double  pdf_area;     /* approximate area below the pdf                    */
  double  mode;         /* location of mode                                  */
  double  il;           /* left border of the domain                         */
  double  ir;           /* right border of the domain                        */
  double *pdf_param;    /* parameters of the pdf                             */
  int     n_pdf_param;  /* number of parameters of the pdf                   */
  double  vollc,volcompl,voll,fm,hm,
    al,ar,col,cor,sal,sar,bl,br,tlx,trx,
    brblvolc,drar,dlal,ooar2,ooal2;/*constants of the hat and for generation*/
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_utdr_new( double (*pdf)(double x,double *pdf_param, int n_pdf_params), double mode );
/* get default parameters for generator                                      */

struct unur_gen *unur_utdr_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_utdr_sample( struct unur_gen *generator );
double unur_utdr_sample_check( struct unur_gen *generator );  /** TODO **/
/* sample from generator                                                     */

void unur_utdr_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*---------------------------------------------------------------------------*/

