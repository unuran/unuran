/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for parameter and generator objects.          *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unur_struct.h                                    *
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
/* types for sampling routines                                               */

/* for univariate continuous distribution */
typedef double _UNUR_SAMPLING_ROUTINE_CONT(struct unur_gen *gen);

/* for univariate discrete distribution */
typedef int _UNUR_SAMPLING_ROUTINE_DISCR(struct unur_gen *gen);

/* for multivariate continuous distribution */
typedef void _UNUR_SAMPLING_ROUTINE_VEC(struct unur_gen *gen, double *vec);


/*---------------------------------------------------------------------------*/
/* parameter objects                                                         */

struct unur_par {
  void *datap;                /* pointer to data for method                  */
  size_t s_datap;             /* size of data structure                      */

  struct unur_gen* (*init)(struct unur_par *par);

  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  unsigned set;               /* stores which parameters have been changed   */

  UNUR_URNG *urng;            /* pointer to uniform random number generator  */
  UNUR_URNG *urng_aux;        /* pointer to second (auxiliary) uniform RNG   */

  const struct unur_distr *distr;   /* pointer to distribution object        */

  unsigned debug;             /* debugging flags                             */
#ifdef UNUR_COOKIES
  unsigned cookie;            /* magic cookie                                */
#endif
};


/*---------------------------------------------------------------------------*/
/* generator objects                                                         */

struct unur_gen { 
  void *datap;                /* pointer to data for method                  */
  
  union {
    _UNUR_SAMPLING_ROUTINE_CONT  *cont;
    _UNUR_SAMPLING_ROUTINE_DISCR *discr;
    _UNUR_SAMPLING_ROUTINE_VEC   *cvec;
    _UNUR_SAMPLING_ROUTINE_VEC   *matr;
  } sample;                   /* pointer to sampling routine                 */
  
  UNUR_URNG *urng;            /* pointer to uniform random number generator  */
  UNUR_URNG *urng_aux;        /* pointer to second (auxiliary) uniform RNG   */

  struct unur_distr *distr;   /* distribution object                         */

  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  unsigned set;               /* stores which parameters have been changed   */
  unsigned status;            /* status of generator object                  */
  
  char *genid;                /* identifier for generator                    */

  struct unur_gen *gen_aux;   /* pointer to auxiliary generator object       */
  struct unur_gen **gen_aux_list; /* list of pointers to auxiliary generator objects */

  size_t s_datap;             /* size of data structure                      */
  unsigned debug;             /* debugging flags                             */

  void (*destroy)(struct unur_gen *gen); /* pointer to destructor            */ 
  struct unur_gen* (*clone)(const struct unur_gen *gen ); /* clone generator */

#ifdef UNUR_COOKIES
  unsigned cookie;            /* magic cookie                                */
#endif
};

/*---------------------------------------------------------------------------*/
