/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_tdr.h                                                        *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method TDR        *
 *         (Transformed Density Rejection)                                   *
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

struct unur_tdr_par { 

  double (*pdf)(double x, double *pdf_param, int n_pdf_param);  /* pointer to pdf */
  double (*dpdf)(double x, double *pdf_param, int n_pdf_param);  /* derivative of pdf */
  double *pdf_param;            /* parameters of the pdf                     */
  int     n_pdf_param;          /* number of parameters of the pdf           */
  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */
  double  mode;                 /* (approximate) location of mode            */

  double  guide_factor;         /* relative size of guide table              */

  double *starting_cpoints;     /* pointer to array of starting points       */
  int     n_starting_cpoints;   /* number of construction points at start    */
  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* bound for ratio r_n = Atotal / Asqueeze   */
  double  bound_for_adding;     /* lower bound for relative area             */

  double  c_T;                  /* parameter c for transformation T_c        */           
};

/*---------------------------------------------------------------------------*/
/* store data for segments                                                   */

struct unur_tdr_interval {

  double  x;                    /* left construction point                   */
  double  fx;                   /* value of p.d.f. at tp                     */ 
  double  Tfx;                  /* value of transformed p.d.f. at tp         */ 
  double  dTfx;                 /* derivative of transformed p.d.f. at  tp   */
  double  sq;                   /* slope of transformed squeeze in interval  */

  double  Acum;                 /* cumulated area of intervals               */
  double  Ahatl;                /* area between hat and squeeze on left side */
  double  Ahatr;                /* area between hat and squeeze on right side*/
  double  Asqueeze;             /* area squeeze                              */

  struct unur_tdr_interval *next; /* pointer to next segment in list         */

#if UNUR_DEBUG & UNUR_DB_COOKIES
  unsigned long cookie;         /* magic cookie                              */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_tdr_gen { 

  double  Atotal;               /* area below hat                            */
  double  Asqueeze;             /* area below squeeze                        */

  double (*pdf)(double x, double *pdf_param, int n_pdf_param);  /* pointer to pdf */
  double (*dpdf)(double x, double *pdf_param, int n_pdf_param);  /* derivative of pdf */
  double *pdf_param;            /* parameters of the pdf                     */
  int     n_pdf_param;          /* number of parameters of the pdf           */
  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */

  double  c_T;                  /* parameter c for transformation T_c        */           

  struct unur_tdr_interval *iv; /* pointer to linked list of intervals       */
  int     n_ivs;                /* number of intervals                       */
  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* bound for ratio r_n = Atotal / Asqueeze   */
  double  bound_for_adding;     /* lower bound for relative area             */

  struct unur_tdr_interval **guide; /* pointer to guide table                */
  int     guide_size;           /* size of guide table                       */
  double  guide_factor;         /* relative size of guide table              */

  struct unur_tdr_interval *iv_stack; /* stack of allocated intervals        */
  int     iv_free;              /* position of last free segment in stack    */
  struct unur_mblock *mblocks;  /* linked list for allocated blocks          */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_tdr_new( double (*pdf)(double x,double *pdf_param, int n_pdf_params), 
			       double (*dpdf)(double x,double *pdf_param, int n_pdf_params) );
/* get default parameters for generator                                      */

struct unur_gen *unur_tdr_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_tdr_sample_log( struct unur_gen *gen );
double unur_tdr_sample_sqrt( struct unur_gen *gen );
double unur_tdr_sample_check( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_tdr_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*---------------------------------------------------------------------------*/

