/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tdr_struct.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method TDR                                *
 *         (Transformed Density Rejection)                                   *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
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

  double  center;               /* approximate location of mode              */

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

#ifdef UNUR_COOKIES
  unsigned cookie;              /* magic cookie                              */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_tdr_gen { 

  double  Atotal;               /* area below hat                            */
  double  Asqueeze;             /* area below squeeze                        */

  double  c_T;                  /* parameter c for transformation T_c        */           

  struct unur_tdr_interval *iv; /* pointer to linked list of intervals       */
  int     n_ivs;                /* number of intervals                       */
  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* bound for ratio r_n = Atotal / Asqueeze   */
  double  bound_for_adding;     /* lower bound for relative area             */

  struct unur_tdr_interval **guide; /* pointer to guide table                */
  int     guide_size;           /* size of guide table                       */
  double  guide_factor;         /* relative size of guide table              */
};

/*---------------------------------------------------------------------------*/
