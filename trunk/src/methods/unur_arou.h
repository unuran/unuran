/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_arou.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method AROU       *
 *         (Adaptive Ratio-Of-Uniforms)                                      *
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

struct unur_arou_par { 

  double (*pdf)(double x, double *pdf_param, int n_pdf_param);  /* pointer to pdf */
  double (*dpdf)(double x, double *pdf_param, int n_pdf_param);  /* derivative of pdf */
  double *pdf_param;            /* parameters of the pdf                     */
  int     n_pdf_param;          /* number of parameters of the pdf           */
  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */
  double  mode;                 /* (approximate) location of mode            */

  double  guide_factor;         /* relative size of guide table              */

  double  bound_for_adding;     /* lower bound for relative area             */
  double  max_ratio;            /* limit for ratio r_n = |P^s| / |P^e|       */
  int     n_starting_cpoints;   /* number of construction points at start    */
  double *starting_cpoints;     /* pointer to array of starting points       */
  int     max_segs;             /* maximum number of segments                */
};

/*---------------------------------------------------------------------------*/
/* store data for segments                                                   */

struct unur_arou_segment {
  double Acum;                  /* cumulated sum of areas                    */
  double Ain;                   /* area of segment inside of squeeze         */
  double Aout;                  /* area of segment outside of squeeze        */

  double ltp[2];                /* coordinates of left tp point in segment   */
  double dltp[3];               /* tanget line of region at left touching point:
				   dltp[0]*u + dltp[1]*v == dltp[2]          */
  double mid[2];                /* coordinates of middle (outer) vertex of segment */
  double *rtp;                  /* pointer to coordinates of right tp in segment
				   (stored in next segment)                  */
  double *drtp;                 /* pointer to tangent line at right tp       */

  struct unur_arou_segment *next; /* pointer to next segment in list         */

#if UNUR_DEBUG & UNUR_DB_COOKIES
  int cookie;                   /* magic cookie                              */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_arou_gen { 

  double  Atotal;               /* area of enveloping polygon                */
  double  Asqueeze;             /* area of squeeze polygon                   */

  double (*pdf)(double x, double *pdf_param, int n_pdf_param);  /* pointer to pdf */
  double (*dpdf)(double x, double *pdf_param, int n_pdf_param);  /* derivative of pdf */
  double *pdf_param;            /* parameters of the pdf                     */
  int     n_pdf_param;          /* number of parameters of the pdf           */
  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */

  double  max_ratio;            /* limit for ratio r_n = |P^s| / |P^e|       */

  struct unur_arou_segment **guide;  /* pointer to guide table               */
  int     guide_size;           /* size of guide table                       */
  double  guide_factor;         /* relative size of guide table              */

  struct unur_arou_segment *seg;     /* pointer to linked list of segments   */
  int     n_segs;               /* number of construction points             */
  int     max_segs;             /* maximum number of segments                */

  double  bound_for_adding;     /* lower bound for relative area             */

  struct unur_arou_segment *seg_stack; /* stack of allocated segments        */
  int     seg_free;             /* position of last free segment in stack    */
  struct unur_mblock  *mblocks; /* linked list for allocated blocks          */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_arou_new( double (*pdf)(double x,double *pdf_param, int n_pdf_params), 
			   double (*dpdf)(double x,double *pdf_param, int n_pdf_params) );
/* get default parameters for generator                                      */

struct unur_gen *unur_arou_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_arou_sample( struct unur_gen *generator );
double unur_arou_sample_check( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_arou_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*---------------------------------------------------------------------------*/
