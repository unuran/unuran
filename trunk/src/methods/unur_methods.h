/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_methods.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and declares structures and function prototypes    *
 *         for all UNURAN methods                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         included in all methods source files.                             *
 *         required for every application of UNURAN.                         *
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
#ifndef __UNUR_METHODS_H_SEEN
#define __UNUR_METHODS_H_SEEN
/*---------------------------------------------------------------------------*/

/** TODO **/
#include <unur_methods_lib.h>

#include <float.h>
#include <stdlib.h>

#include <unur_defs.h>

#include <unur_umalloc.h>
#include <unur_urng.h>

/*---------------------------------------------------------------------------*/
/* Typedefs                                                                  */

struct unur_gen;
struct unur_par;

typedef double _UNUR_FUNCTION_CONT(double x, double *params, int n_params);
/* a function for continuous univariate c.d.f., p.d.f. and its derivative    */

/* sampling routines                                                         */
typedef double _UNUR_SAMPLING_ROUTINE_CONT(struct unur_gen *gen);
/* for univariate continuous distribution                                    */

typedef int _UNUR_SAMPLING_ROUTINE_DISCR(struct unur_gen *gen);
/* for univariate discrete distribution                                      */

typedef void _UNUR_SAMPLING_ROUTINE_VEC(struct unur_gen *gen, double *vec);
/* for multivariate continuous distribution                                  */

/*---------------------------------------------------------------------------*/
/* include header file for distribtion object                                */
#include <unur_distr.h>

/*---------------------------------------------------------------------------*/
/* include header files for generators                                       */

/* methods for discrete distributions                                        */
#include <unur_dau.h>
#include <unur_dis.h>

/* methods for continuous distributions                                      */
#include <unur_arou.h>
#include <unur_srou.h>
#include <unur_stdr.h>
#include <unur_tabl.h>
#include <unur_tdr.h>
#include <unur_unif.h>
#include <unur_utdr.h>

/* methods for continuous multivariate distributions                         */
#include <unur_rect.h>

/* generators for standard distributions                                     */
#include <unur_cstd.h>

/*---------------------------------------------------------------------------*/
/* get type of transformation method                                         */

#define unur_is_discr(gen) ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ? 1 : 0 )
#define unur_is_cont(gen)  ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_CONT)  ? 1 : 0 )
#define unur_is_vec(gen)   ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_VEC)   ? 1 : 0 )

/*---------------------------------------------------------------------------*/
/* Main structure for all UNURAN generators                                  */

/* parameters */
struct unur_par {
  union {             
    struct unur_dau_par   dau;
    struct unur_dis_par   dis;
    struct unur_arou_par  arou;
    struct unur_srou_par  srou;
    struct unur_stdr_par  stdr;
    struct unur_tabl_par  tabl;
    struct unur_tdr_par   tdr;
    struct unur_unif_par  unif;
    struct unur_utdr_par  utdr;
    struct unur_rect_par  rect;
    struct unur_cstd_par  cstd;
  }               data;       /* data for method                             */

  struct unur_gen* (*init)(struct unur_par *par);

  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  unsigned set;               /* stores which parameters have been changed   */

  UNUR_URNG_TYPE  urng;       /* pointer to uniform random number generator  */

  struct unur_distr *distr;   /* pointer to distribution object              */

  unsigned debug;             /* debugging flags                             */
#if UNUR_DEBUG & UNUR_DB_COOKIES  /* use magic cookies */
  unsigned cookie;            /* magic cookie                                */
#endif
};

/* generators */
struct unur_gen { 
  union {   
    struct unur_dau_gen   dau;
    struct unur_dis_gen   dis;
    struct unur_arou_gen  arou;
    struct unur_srou_gen  srou;
    struct unur_stdr_gen  stdr;
    struct unur_tabl_gen  tabl;
    struct unur_tdr_gen   tdr;
    struct unur_unif_gen  unif;
    struct unur_utdr_gen  utdr;
    struct unur_rect_gen  rect;
    struct unur_cstd_gen  cstd;
  }               data;       /* data for method                             */
  
  union {
    _UNUR_SAMPLING_ROUTINE_CONT  *cont;
    _UNUR_SAMPLING_ROUTINE_DISCR *discr;
    _UNUR_SAMPLING_ROUTINE_VEC   *vec;
  }               sample;     /* pointer to sampling routine                 */
  
  void (*destroy)(struct unur_gen *gen); /* pointer to destructor            */ 
  
  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  
  UNUR_URNG_TYPE  urng;       /* pointer to uniform random number generator  */

  struct unur_distr distr;    /* distribution object                         */
  char *genid;                /* identifier for generator                    */
  
  unsigned debug;             /* debugging flags                             */
#if UNUR_DEBUG & UNUR_DB_COOKIES  /* use magic cookies */
  unsigned cookie;            /* magic cookie                                */
#endif
};

/*---------------------------------------------------------------------------*/
/* invoke generators                                                         */

#if UNUR_DEBUG & UNUR_DB_CHECKNULL
/* check for (invalid) NULL pointer */

#define unur_init(par)                ((par) ? (par)->init(par) : NULL)

#define unur_sample_discr(gen)        ((gen) ? (gen)->sample.discr(gen) : 0)
#define unur_sample_cont(gen)         ((gen) ? (gen)->sample.cont(gen)  : 0.)
#define unur_sample_vec(gen,vector)   if(gen) (gen)->sample.vec(gen,vector)

#define unur_free(gen)                if (gen) (gen)->destroy(gen)

#else

#define unur_init(par)                par->init(par)

#define unur_sample_discr(gen)        (gen)->sample.discr(gen)
#define unur_sample_cont(gen)         (gen)->sample.cont(gen)
#define unur_sample_vec(gen,vector)   (gen)->sample.vec(gen,vector)

#define unur_free(gen)                gen->destroy(gen)

#endif  /* UNUR_DB_CHECKNULL */

/*---------------------------------------------------------------------------*/
/* misc prototype                                                            */
#include <unur_misc.h>

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_METHODS_H_SEEN */
/*---------------------------------------------------------------------------*/


