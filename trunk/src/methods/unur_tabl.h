/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_tabl.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method TABL       *
 *         (Ahren's TABLe method: piecewise constant hat)                    *
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

struct unur_tabl_par { 

  double *slopes;       /* slopes <a_i,b_i>, i.e.\ f(a_i) >= f(b_i)          */
  int     n_slopes;     /* number of slopes                                  */
  double  bleft;        /* left border of the domain                         */
  double  bright;       /* right border of the domain                        */

  int     max_ivs;      /* maximum number of intervals                       */
  double  max_ratio;    /* limit for ratio r_n = A(squeeze) / A(hat)         */

  int     n_starting_cpoints;   /* number of construction points at start    */
  double  area_fract;   /* parameter for equal area rule                     */

  double  guide_factor; /* relative size of guide table                      */
};

/*---------------------------------------------------------------------------*/
/* storing information about generator                                       */

struct unur_tabl_interval {

  double  xmax;         /* maximum of pdf in interval                        */
  double  fmax;         /* maximal value of pdf in interval                  */
  double  xmin;         /* minimum of pdf in interval                        */
  double  fmin;         /* minimal value of pdf in interval                  */
  int     slope;        /* decreasing (-1) or increasing (+1)    (= sign(f´(x))
			   <a,b> = [a,b]      <a,b> = [b,a]                  */

  double  Ahat;         /* area of total bar (below hat)                     */
  double  Asqueeze;     /* area of bar below squeeze                         */
  double  Acum;         /* cumulated area of bars                            */

  struct unur_tabl_interval *next;  /* pointer to next element in list       */

#if UNUR_DEBUG & UNUR_DB_COOKIES
  unsigned cookie;      /* magic cookie                                      */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_tabl_gen { 
  double  Atotal;               /* total area below hat                      */
  double  Asqueeze;             /* area of squeeze polygon                   */

  _UNUR_FUNCTION_CONT *pdf;     /* pointer to p.d.f.                         */
  double *pdf_param;            /* parameters of the pdf                     */
  int     n_pdf_param;          /* number of parameters of the pdf           */
  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */

  struct unur_tabl_interval **guide; /* pointer to guide table               */
  int     guide_size;           /* size of guide table                       */
  double  guide_factor;         /* relative size of guide table              */

  struct unur_tabl_interval *iv;     /* pointer to linked list of intervals  */
  int     n_ivs;                /* number of intervals                       */
  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* limit for ratio r_n = A(squeeze) / A(hat) */

  struct unur_tabl_interval *iv_stack; /* stack of allocated intervals       */
  int     iv_free;              /* position of last free segment in stack    */
  struct unur_mblock *mblocks;  /* linked list for allocated blocks          */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_tabl_new( struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_tabl_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_tabl_sample( struct unur_gen *generator );
double unur_tabl_sample_check( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_tabl_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_tabl_set_nstp( struct unur_par *par, int n_stp );
/* set number of construction points for hat at initialization               */

int unur_tabl_set_max_sqhratio( struct unur_par *par, double max_ratio );
/* set bound for ratio A(squeeze) / A(hat)                                   */

int unur_tabl_set_max_intervals( struct unur_par *par, int max_ivs );
/* set maximum number of intervals                                           */

int unur_tabl_set_areafraction( struct unur_par *par, double fraction );
/* set parameter for equal area rule                                         */

int unur_tabl_set_slopes( struct unur_par *par, double *slopes, int n_slopes );
/* set slopes of p.d.f.                                                      */

int unur_tabl_set_boundary( struct unur_par *par, double left, double right );
/* set left and right boundary of computation interval                       */

int unur_tabl_set_variant( struct unur_par *parameters, unsigned variant );
/* set variant for method                                                    */

int unur_tabl_set_guidefactor( struct unur_par *par, double factor );
/* set factor for relative size of guide table                               */

int unur_tabl_set_verify( struct unur_par *par, int verify );
/* turn verifying of algorithm while sampling on/off                         */

#define unur_dis_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/



