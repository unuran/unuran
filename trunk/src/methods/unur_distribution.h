/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_distribution.h                                               *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines distribution object and                                   *
 *         declares function prototypes for manipulating such an object      *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only be used in unur_methods.h and unur_distr.h                   *
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
#ifndef __UNUR_DISTRIBUTION_H_SEEN
#define __UNUR_DISTRIBUTION_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Typedefs                                                                  */

/* function that return pointer to sampling routine                          */
typedef _UNUR_SAMPLING_ROUTINE_CONT *_UNUR_GET_SAMPLING_ROUTINE_CONT(unsigned);
/* for univariate continuous distribution                                    */

/* function that returns name of sampling routine                            */
typedef const char *_UNUR_GET_SAMPLING_NAME(void *);

/*---------------------------------------------------------------------------*/
/* define object for univariate continuous distribution                      */
struct unur_distr_cont {

  _UNUR_FUNCTION_CONT *pdf;     /* pointer to p.d.f.                         */
  _UNUR_FUNCTION_CONT *dpdf;    /* pointer to derivative of p.d.f.           */
  _UNUR_FUNCTION_CONT *cdf;     /* pointer to c.d.f.                         */

  double params[UNUR_DISTR_MAXPARAMS];  /* parameters of the p.d.f.          */
  int    n_params;              /* number of parameters of the pdf           */

  double mode;                  /* location of mode                          */
  double area;                  /* area below p.d.f.                         */
  double domain[2];             /* boundary of domain                        */

  _UNUR_GET_SAMPLING_ROUTINE_CONT *get_sampling_routine;  /* get pointer to sampling routine */
#if UNUR_DEBUG & UNUR_DB_INFO
  _UNUR_GET_SAMPLING_NAME         *get_sampling_name;     /* get name of sampling routine */
#endif
};

/*---------------------------------------------------------------------------*/
/* define object for univariate discrete distribution                        */
struct unur_distr_discr {

  double *prob;                 /* pointer to probability vector             */
  int     n_prob;               /* length of probability vector              */
};

/*---------------------------------------------------------------------------*/
/* define distribution object                                                */

struct unur_distr {
  union {             
    struct unur_distr_cont  cont;   /* univariate continuous distribution    */
    struct unur_distr_discr discr;  /* univariate discrete distribution      */
  } data;                           /* data for distribution                 */

  unsigned type;                    /* type of distribution                  */
  unsigned id;                      /* identifier for distribution           */
  const char *name;                 /* name of distribution                  */

  unsigned set;                     /* indicate changed parameters           */

#if UNUR_DEBUG & UNUR_DB_COOKIES    /* use magic cookies */
  unsigned cookie;                  /* magic cookie                          */
#endif
};

/*---------------------------------------------------------------------------*/
/* types of distribtuions                                                    */

enum {
  UNUR_DISTR_CONT  = 0x001u,        /* univariate continuous distribution    */ 
  UNUR_DISTR_DISCR = 0x002u,        /* univariate discrete distribution      */ 
};

/*---------------------------------------------------------------------------*/
/* indicate changed parameters                                               */

enum {
  UNUR_DISTR_SET_PARAMS     = 0x001u,
  UNUR_DISTR_SET_DOMAIN     = 0x002u,
  UNUR_DISTR_SET_STDDOMAIN  = 0x004u,   /* domain not truncated (for standard distributions) */
  UNUR_DISTR_SET_MODE       = 0x008u,
  UNUR_DISTR_SET_PDFAREA    = 0x010u,
}; 

/*---------------------------------------------------------------------------*/
/* indentifiers for standard distributions                                   */

enum {
  UNUR_DISTR_GENERIC  = 0x0u,

  UNUR_DISTR_BETA,
  UNUR_DISTR_CAUCHY,
  UNUR_DISTR_CHISQUARE,
  UNUR_DISTR_EXPONENTIAL,
  UNUR_DISTR_GAMMA,
  UNUR_DISTR_LAPLACE,
  UNUR_DISTR_LOGNORMAL,
  UNUR_DISTR_LOMAX,
  UNUR_DISTR_NORMAL,
  UNUR_DISTR_PARETO,
  UNUR_DISTR_POWEREXPONENTIAL,
  UNUR_DISTR_RAYLEIGH,
  UNUR_DISTR_STUDENT,
  UNUR_DISTR_UNIFORM,
};

/*---------------------------------------------------------------------------*/
/* Common variants for all special generators                                */

#define UNUR_STDGEN_INVERSION (~0u)

/*---------------------------------------------------------------------------*/
/* function prototypes for manipulating distribution object                  */

/*---------------------------------------------------------------------------*/
/* routines for all distribution objects                                     */

struct unur_distr *unur_distr_dup( struct unur_distr *distr );
/* duplicate distribution object                                             */

void unur_distr_free( struct unur_distr *distr );
/* destroy distribution object                                               */

int unur_distr_set_name( struct unur_distr *distr, const char *name );
/* set name of distribution                                                  */

/*---------------------------------------------------------------------------*/
/* univariate continuous distributions                                       */

struct unur_distr *unur_distr_cont_new( void );
/* create a new distribution object                                          */

int unur_distr_cont_set_pdf( struct unur_distr *distr, void *pdf );
/* set p.d.f. of distribution                                                */

int unur_distr_cont_set_dpdf( struct unur_distr *distr, void *dpdf );
/* set derivative of p.d.f. of distribution                                  */

int unur_distr_cont_set_cdf( struct unur_distr *distr, void *cdf );
/* set c.d.f. of distribution                                                */

int unur_distr_cont_set_params( struct unur_distr *distr, double *params, int n_params );
/* set array of parameters for distribution                                  */

int unur_distr_cont_set_mode( struct unur_distr *distr, double mode );
/* set mode of distribution                                                  */

int unur_distr_cont_set_pdfarea( struct unur_distr *distr, double area );
/* set area below p.d.f.                                                     */

int unur_distr_cont_set_domain( struct unur_distr *distr, double left, double right );
/* set the left and right borders of the domain of the distribution          */

void _unur_distr_cont_debug( struct unur_distr *distr, char *genid );
/* write info about distribution into logfile                                */

/*---------------------------------------------------------------------------*/
/* discrete univariate distributions                                         */

struct unur_distr *unur_distr_discr_new( void );
/* create a new distribution object                                          */

int unur_distr_discr_set_prob( struct unur_distr *distr, double *prob, int n_prob );
/* set probability vector for distribution                                   */

void _unur_distr_discr_debug( struct unur_distr *distr, char *genid, int printvector );
/* write info about distribution into logfile                                */

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_DISTRIBUTION_H_SEEN */
/*---------------------------------------------------------------------------*/

