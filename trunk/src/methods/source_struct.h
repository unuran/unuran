/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_struct.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for distribution, parameter, and generator    *
 *         objects.                                                          *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_unuran.h                                  *
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
#ifndef SOURCE_STRUCT_H_SEEN
#define SOURCE_STRUCT_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Basic header files                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* compiler switches and defaults */
#include <unuran_config.h>

/*****************************************************************************/
/**  Functions                                                              **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* the objects                                                               */
struct unur_distr;    /* distribution object      */
struct unur_par;      /* parameters for generator */
struct unur_gen;      /* generator object         */

/*---------------------------------------------------------------------------*/
/* sampling routines                                                         */

/* for univariate continuous distribution */
typedef double _UNUR_SAMPLING_ROUTINE_CONT(struct unur_gen *gen);

/* for univariate discrete distribution */
typedef int _UNUR_SAMPLING_ROUTINE_DISCR(struct unur_gen *gen);

/* for multivariate continuous distribution */
typedef void _UNUR_SAMPLING_ROUTINE_VEC(struct unur_gen *gen, double *vec);

/*****************************************************************************/
/**  Auxiliary tools                                                        **/
/*****************************************************************************/

#include <x_slist_struct.h>

/*****************************************************************************/
/**  Declaration for parser                                                 **/
/*****************************************************************************/

#include <functparser_struct.h>

/*****************************************************************************/
/**  Declarations for uniform random number generators                      **/
/*****************************************************************************/

#include <x_urng.h>

/*****************************************************************************/
/**  Distribution objects                                                   **/
/*****************************************************************************/

#include <distr_struct.h>

/*****************************************************************************/
/**  structures for generators                                              **/
/*****************************************************************************/

/* automatically selected method */
#include <auto_struct.h>

/* discrete distributions */
#include <dari_struct.h>
#include <dau_struct.h>
#include <dgt_struct.h>
#include <dsrou_struct.h>

/* continuous distributions */
#include <arou_struct.h>
#include <ninv_struct.h>
#include <srou_struct.h>
#include <ssr_struct.h>
#include <tabl_struct.h>
#include <tdr_struct.h>
#include <unif_struct.h>
#include <utdr_struct.h>

#include <empk_struct.h>

/* continuous multivariate distributions */
#include <vempk_struct.h>
#include <vmt_struct.h>

/* wrappers for special generators for standard distributions */
#include <cstd_struct.h>     /* continuous */
#include <dstd_struct.h>     /* discrete   */

/*****************************************************************************/
/**  Main structure for all UNURAN generators                               **/  
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* parameter objects                                                         */

struct unur_par {
  union {             
    struct unur_auto_par  mauto;
    struct unur_dari_par  dari;
    struct unur_dau_par   dau;
    struct unur_dgt_par   dgt;
    struct unur_dsrou_par dsrou;
    struct unur_arou_par  arou;
    struct unur_ninv_par  ninv;
    struct unur_srou_par  srou;
    struct unur_ssr_par   ssr;
    struct unur_tabl_par  tabl;
    struct unur_tdr_par   tdr;
    struct unur_unif_par  unif;
    struct unur_utdr_par  utdr;
    struct unur_empk_par  empk;
    struct unur_vmt_par   vmt;
    struct unur_vempk_par vempk;
    struct unur_cstd_par  cstd;
    struct unur_dstd_par  dstd;
  }               data;       /* data for method                             */

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
  union {   
    struct unur_auto_gen  mauto;
    struct unur_dari_gen  dari;
    struct unur_dau_gen   dau;
    struct unur_dgt_gen   dgt;
    struct unur_dsrou_gen dsrou;
    struct unur_arou_gen  arou;
    struct unur_ninv_gen  ninv;
    struct unur_srou_gen  srou;
    struct unur_ssr_gen   ssr;
    struct unur_tabl_gen  tabl;
    struct unur_tdr_gen   tdr;
    struct unur_unif_gen  unif;
    struct unur_utdr_gen  utdr;
    struct unur_empk_gen  empk;
    struct unur_vmt_gen   vmt;
    struct unur_vempk_gen vempk;
    struct unur_cstd_gen  cstd;
    struct unur_dstd_gen  dstd;
  }               data;       /* data for method                             */
  
  union {
    _UNUR_SAMPLING_ROUTINE_CONT  *cont;
    _UNUR_SAMPLING_ROUTINE_DISCR *discr;
    _UNUR_SAMPLING_ROUTINE_VEC   *cvec;
  }               sample;     /* pointer to sampling routine                 */
  
  void (*destroy)(struct unur_gen *gen); /* pointer to destructor            */ 
  
  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  unsigned set;               /* stores which parameters have been changed   */
  
  UNUR_URNG *urng;            /* pointer to uniform random number generator  */
  UNUR_URNG *urng_aux;        /* pointer to second (auxiliary) uniform RNG   */

  struct unur_distr distr;    /* distribution object                         */
  char *genid;                /* identifier for generator                    */

  struct unur_gen *gen_aux;   /* pointer to auxiliary generator object       */

  unsigned debug;             /* debugging flags                             */

#ifdef UNUR_COOKIES
  unsigned cookie;            /* magic cookie                                */
#endif
};

/*---------------------------------------------------------------------------*/
#endif  /* SOURCE_STRUCT_H_SEEN */
/*---------------------------------------------------------------------------*/
