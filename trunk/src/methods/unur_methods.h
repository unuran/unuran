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

#include <float.h>
#include <stdlib.h>

#include <unur_defs.h>

#include <unur_umalloc.h>
#include <unur_urng.h>

/*---------------------------------------------------------------------------*/
/* include header file for distribtion object                                */
#include <unur_distribution.h>

/*---------------------------------------------------------------------------*/
/* include header files for generators                                       */

/* methods for discrete distributions                                        */
#include <unur_dau.h>
#include <unur_dis.h>

/* methods for continuous distributions                                      */
#include <unur_arou.h>
#include <unur_tabl.h>
#include <unur_tdr.h>
#include <unur_unif.h>
#include <unur_utdr.h>

/* methods for continuous multivariate distributions                         */
#include <unur_rect.h>

/* generators for standard distributions                                     */
#include <unur_cstd.h>

/* set parameters for generators                                             */
#include <unur_set.h>

/*---------------------------------------------------------------------------*/
/* List of methods                                                           */

#define UNUR_MASK_VARIANT  0x00000fffUL   /* indicate variant (see the corresponding .c files) */
#define UNUR_MASK_METHOD   0xfff00000UL   /* indicate method                   */
#define UNUR_MASK_TYPE     0xf0000000UL   /* indicate type of method           */

/* bits 13-20 are used for flags common to all generators */
#define UNUR_MASK_SCHECK   0x00001000UL   /* turns check sampling on/off       */
#define UNUR_MASK_MODE     0x00002000UL   /* use mode                          */
#define UNUR_MASK_COPYALL  0x00004000UL   /* copy all inputs into generator object */

/* discrete, univariate */
#define UNUR_METH_DISCR    0x10000000UL

#define UNUR_METH_DAU      0x10100000UL
#define UNUR_METH_DIS      0x10200000UL

/* continuous, univariate */
#define UNUR_METH_CONT     0x20000000UL

#define UNUR_METH_AROU     0x20300000UL
#define UNUR_METH_TABL     0x20400000UL
#define UNUR_METH_TDR      0x20500000UL
#define UNUR_METH_UNIF     0x20600000UL
#define UNUR_METH_UTDR     0x20700000UL

/* continuous, multivariate */
#define UNUR_METH_VEC      0x40000000UL

#define UNUR_METH_RECT     0x40700000UL

/* generators for standard distributions                                     */
/* for definitions of methods for standard distributions see "stand.c"       */
#define UNUR_MASK_DISTR    0x000ffff0UL   /* indicate distribution           */

#define UNUR_METH_CSTD     0x2f000000UL   /* is of type UNUR_METH_CONT !! */

/* to indicate unkown type */
#define UNUR_METH_UNKNOWN  0xf0000000UL


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
    struct unur_tabl_par  tabl;
    struct unur_tdr_par   tdr;
    struct unur_unif_par  unif;
    struct unur_utdr_par  utdr;
    struct unur_rect_par  rect;
    struct unur_cstd_par  cstd;
  }               data;       /* data for method                             */

  struct unur_gen* (*init)(struct unur_par *par);

  unsigned long   method;     /* indicates method and generator to be used   */
  unsigned long   set;        /* stores which parameters have been changed   */

  UNUR_URNG_TYPE  urng;       /* pointer to uniform random number generator  */

  struct unur_distr *distr;   /* pointer to distribution object              */

#if UNUR_DEBUG & UNUR_DB_COOKIES  /* use magic cookies */
  unsigned long   cookie;     /* magic cookie                                */
#endif
#if UNUR_DEBUG & UNUR_DB_INFO     /* print data about generators */
  unsigned long   debug;      /* debugging flags                             */
#endif
};

/* generators */
struct unur_gen { 
  union {   
    struct unur_dau_gen   dau;
    struct unur_dis_gen   dis;
    struct unur_arou_gen  arou;
    struct unur_tabl_gen  tabl;
    struct unur_tdr_gen   tdr;
    struct unur_unif_gen  unif;
    struct unur_utdr_gen  utdr;
    struct unur_rect_gen  rect;
    struct unur_cstd_gen  cstd;
  }               data;       /* data for method                             */
  
  union {
    int    (*discr)(struct unur_gen *gen);
    double (*cont) (struct unur_gen *gen);
    void   (*vec)  (struct unur_gen *gen, double *vec);
  }               sample;     /* pointer to sampling routine                 */
  
  void            (*destroy)(struct unur_gen *gen); /* pointer to destructor */ 
  
  unsigned long   method;     /* indicates method and generator to be used   */
  
  UNUR_URNG_TYPE  urng;       /* pointer to uniform random number generator  */

  struct unur_distr *distr;   /* pointer to distribution object              */
  
#if UNUR_DEBUG & UNUR_DB_COOKIES  /* use magic cookies */
  unsigned long   cookie;     /* magic cookie                                */
#endif
#if UNUR_DEBUG & UNUR_DB_INFO     /* print data about generators */
  unsigned long   debug;      /* debugging flags                             */
#endif
#if UNUR_DEBUG > 0                /* debugging compiled into generators */
  char           *genid;      /* identifier for generator                    */
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
#endif  /* __UNUR_METHODS_H_SEEN */
/*---------------------------------------------------------------------------*/

