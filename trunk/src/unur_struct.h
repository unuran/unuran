/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for distribution, parameter, and generator    *
 *         objects.                                                          *
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
#ifndef UNUR_STRUCT_H_SEEN
#define UNUR_STRUCT_H_SEEN
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

/*---------------------------------------------------------------------------*/
/* Generic functions                                                          */

typedef double UNUR_FUNCT_GENERIC  (double  x, void *params);
typedef double UNUR_FUNCT_VGENERIC (double *x, void *params);

/* for univariate functions with optional parameter array */
struct unur_funct_generic {
  UNUR_FUNCT_GENERIC *f;
  void *params;
};

/* for multivariate functions with optional parameter array */
struct unur_funct_vgeneric {
  UNUR_FUNCT_VGENERIC *f;
  void *params;
};

/*****************************************************************************/
/**  Auxiliary tools                                                        **/
/*****************************************************************************/

#include <utils/slist_struct.h>
#include <utils/string_struct.h>

/*****************************************************************************/
/**  Declaration for parser                                                 **/
/*****************************************************************************/

#include <parser/functparser_struct.h>

/*****************************************************************************/
/**  Declarations for uniform random number generators                      **/
/*****************************************************************************/

#include <methods/x_urng.h>

/*****************************************************************************/
/**  Distribution objects                                                   **/
/*****************************************************************************/

#include <distr/cvec_struct.h>
#include <distr/distr_struct.h>

/*****************************************************************************/
/**  structures for generators                                              **/
/*****************************************************************************/

/* automatically selected method */
#include <methods/auto_struct.h>

/* discrete distributions */
#include <methods/dari_struct.h>
#include <methods/dau_struct.h>
#include <methods/dgt_struct.h>
#include <methods/dsrou_struct.h>
#include <methods/dss_struct.h>

/* continuous distributions */
#include <methods/arou_struct.h>
#include <methods/hinv_struct.h>
#include <methods/hrb_struct.h>
#include <methods/hrd_struct.h>
#include <methods/hri_struct.h>
#include <methods/ninv_struct.h>
#include <methods/nrou_struct.h>
#include <methods/srou_struct.h>
#include <methods/ssr_struct.h>
#include <methods/tabl_struct.h>
#include <methods/tdr_struct.h>
#include <methods/unif_struct.h>
#include <methods/utdr_struct.h>

#include <methods/empk_struct.h>
#include <methods/empl_struct.h>

/* continuous multivariate distributions */
#include <methods/vempk_struct.h>
#include <methods/vmt_struct.h>

/* random matrices */
#include <methods/mcorr_struct.h>

/* wrappers for special generators for standard distributions */
#include <methods/cstd_struct.h>     /* continuous */
#include <methods/dstd_struct.h>     /* discrete   */

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
    struct unur_dss_par   dss;
    struct unur_arou_par  arou;
    struct unur_hinv_par  hinv;
    struct unur_hrb_par   hrb;
    struct unur_hrd_par   hrd;
    struct unur_hri_par   hri;
    struct unur_mcorr_par mcorr;
    struct unur_ninv_par  ninv;
    struct unur_nrou_par  nrou;
    struct unur_srou_par  srou;
    struct unur_ssr_par   ssr;
    struct unur_tabl_par  tabl;
    struct unur_tdr_par   tdr;
    struct unur_unif_par  unif;
    struct unur_utdr_par  utdr;
    struct unur_empk_par  empk;
    struct unur_empl_par  empl;
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
    struct unur_dss_gen   dss;
    struct unur_arou_gen  arou;
    struct unur_hinv_gen  hinv;
    struct unur_hrb_gen   hrb;
    struct unur_hrd_gen   hrd;
    struct unur_hri_gen   hri;
    struct unur_mcorr_gen mcorr;
    struct unur_ninv_gen  ninv;
    struct unur_nrou_gen  nrou;
    struct unur_srou_gen  srou;
    struct unur_ssr_gen   ssr;
    struct unur_tabl_gen  tabl;
    struct unur_tdr_gen   tdr;
    struct unur_unif_gen  unif;
    struct unur_utdr_gen  utdr;
    struct unur_empk_gen  empk;
    struct unur_empl_gen  empl;
    struct unur_vmt_gen   vmt;
    struct unur_vempk_gen vempk;
    struct unur_cstd_gen  cstd;
    struct unur_dstd_gen  dstd;
  }               data;       /* data for method                             */
  
  union {
    _UNUR_SAMPLING_ROUTINE_CONT  *cont;
    _UNUR_SAMPLING_ROUTINE_DISCR *discr;
    _UNUR_SAMPLING_ROUTINE_VEC   *cvec;
    _UNUR_SAMPLING_ROUTINE_VEC   *matr;
  }               sample;     /* pointer to sampling routine                 */
  
  unsigned method;            /* indicates method and generator to be used   */
  unsigned variant;           /* indicates variant of method                 */
  unsigned set;               /* stores which parameters have been changed   */
  
  UNUR_URNG *urng;            /* pointer to uniform random number generator  */
  UNUR_URNG *urng_aux;        /* pointer to second (auxiliary) uniform RNG   */

  struct unur_distr *distr;   /* distribution object                         */
  char *genid;                /* identifier for generator                    */

  struct unur_gen *gen_aux;   /* pointer to auxiliary generator object       */
  struct unur_gen **gen_aux_list; /* list of pointers to auxiliary generator objects */

  unsigned debug;             /* debugging flags                             */

  void (*destroy)(struct unur_gen *gen); /* pointer to destructor            */ 
  struct unur_gen* (*clone)(const struct unur_gen *gen ); /* clone generator */

#ifdef UNUR_COOKIES
  unsigned cookie;            /* magic cookie                                */
#endif
};

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_STRUCT_H_SEEN */
/*---------------------------------------------------------------------------*/
