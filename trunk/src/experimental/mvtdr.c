/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvtdr.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    multivariate transformed density rejection                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given (logarithm of the) PDF of a log-concave distribution;          *
 *      produce a value x consistent with its density.                       *
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
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Leydold, J. (1998): A Rejection Technique for Sampling from         *
 *       Log-Concave Multivariate Distributions,                             *
 *       ACM TOMACS 8(3), pp. 254-280.                                       *
 *                                                                           *
 *   [2] Hoermann, W., J. Leydold, and G. Derflinger (2004):                 *
 *       Automatic Nonuniform Random Variate Generation, Springer, Berlin.   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  VMT generates random vectors for distributions with given mean           *
 *  vector mu and covariance matrix Sigma. It produces random vectors        *
 *  of the form X = L Y + mu, where L is the Cholesky factor of Sigma,       *
 *  i.e. L L^t = Sigma, and Y has independent components of the same         *
 *  distribution with mean 0 and standard deviation 1.                       *
 *                                                                           *
 *  See [2], Sect.11.1.6, Alg.11.3.                                          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distributions/unur_distributions.h>
#include <utils/fmax_source.h>
#include <uniform/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "tdr.h"
#include "mvtdr.h"
#include "mvtdr_struct.h"

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* relative size of guide table compared to number of cones */
#define GUIDE_TABLE_SIZE    1

/* find proper touching point using Brent's algorithm */
#define FIND_TP_START       1.      /* starting point                       */
#define FIND_TP_STEPS_R     10      /* max number of steps --> Infinity     */
#define FIND_TP_STEPS_L     10      /* max number of steps --> 0            */
#define FIND_TP_STEP_SIZE   2.      /* stepsize for searching routine       */
/* find bracket for Brent's algorithm */
#define FIND_TP_STEPS_LB    10      /* max steps for finding lower bound    */
#define FIND_TP_STEPS_UB    10      /* max stepsfor finding upper bound     */
/* acceptable tolerance for Brent's algorithm */
#define FIND_TP_TOL         0.001 



/** TODO **/

/* control triangulation of cones                                  */
#define MAX_N_CONES         10000    /* maximum number of cones (at least 2^(N+T_STEPS_MIN) */
#define T_STEPS_MIN         5        /* minimum number of triangulation steps */
#define OPTIMAL_TP_STEP     100      /* triangulation step when optimal touching points is calculated */

/** fine tuning of generator                                       **/
/* a number is considered to be zero if abs is below this bound     */
#define TOLERANCE   1.e-8

/* move mode to boundary if | mode - boundary | / length < MODE_TO_BOUNDARY */
#define MODE_TO_BOUNDARY    1.E-2

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MVTDR"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static int _unur_mvtdr_simplex_sample( const struct unur_gen *gen, double *U );
/*---------------------------------------------------------------------------*/
/* sample point uniformly on standard simplex.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_gammagen( struct unur_gen *gen, double alpha );
/*---------------------------------------------------------------------------*/
/* create a gamma random variate generator with shape parameter alpha.       */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/* Hat.                                                                      */
/*****************************************************************************/

static int _unur_mvtdr_create_hat( struct unur_gen *gen );
/* compute cones and hat function */


/*****************************************************************************/
/* CONES.                                                                    */
/*****************************************************************************/

static int _unur_mvtdr_initial_cones( struct unur_gen *gen );
/* get initial cones */

static CONE *_unur_mvtdr_cone_new( struct unur_gen *gen );
/* get new (empty) vertex object */

static int _unur_mvtdr_cone_center( struct unur_gen *gen, CONE *c );
/* computer center of cone */

static int _unur_mvtdr_cone_params( struct unur_gen *gen, CONE *c );
/* compute parameters for hat for a cone   (expect touching point and volume below hat) */

static double _unur_mvtdr_cone_logH( struct unur_gen *gen, CONE *c );
/* calculate log of volume below hat for given touching point */

static int _unur_mvtdr_cone_split( struct unur_gen *gen, CONE *c, int step );
/* split a cone */

static int _unur_mvtdr_triangulate(struct unur_gen *gen, int step, int all);
/* make one triangulation step */

/*****************************************************************************/
/* optimal distance for touching points                                      */
/*****************************************************************************/

static double _unur_mvtdr_tp_min( double t, void *p );
/* wrapper for _unur_mvtdr_cone_hatvolume();
   sets cone->tp;
   funtion that must be minimized for optimal touching point */

static int _unur_mvtdr_tp_find( struct unur_gen *gen, CONE *c );
/* find optimal touching point for cone */

static int _unur_mvtdr_tp_search( struct unur_gen *gen, TP_ARG *a );
/* search for proper touching point */

static int _unur_mvtdr_tp_bracket( struct unur_gen *gen, TP_ARG *a );
/* search for proper bracket of minimum of tp_f2min() */

/*****************************************************************************/
/* VERTICES.                                                                 */
/*****************************************************************************/

static int _unur_mvtdr_initial_vertices( struct unur_gen *gen );
/* get vertices of initial cones */

static VERTEX *_unur_mvtdr_vertex_new( struct unur_gen *gen );
/* get new (empty) vertex object */

static int _unur_mvtdr_number_vertices( struct unur_gen *gen, int level );
/* number of vertices in given triangulation level */

static VERTEX *_unur_mvtdr_vertex_on_edge( struct unur_gen *gen, VERTEX **vl );
/* compute new vertex on edge */

/*****************************************************************************/
/* hash table for storing EDGES.                                             */
/*****************************************************************************/

static int _unur_mvtdr_etable_new( struct unur_gen *gen, int size );
/* make new hash table */

static void _unur_mvtdr_etable_free( struct unur_gen *gen );
/* free hash table */

static VERTEX *_unur_mvtdr_etable_find_or_insert( struct unur_gen *gen, VERTEX **vidx );
/* find or insert entry in hash table, return pointer to vertex */

/*****************************************************************************/

static int _unur_mvtdr_make_guide_table( struct unur_gen *gen );
/* create guide table */

/*****************************************************************************/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_mvtdr_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mvtdr_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */
#define dPDF(r,x) _unur_cvec_dPDF((r),(x),(gen->distr))    /* call to dPDF        */

/* an auxiliary generator for gamma variates */
#define GEN_GAMMA  gen->gen_aux

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
  error codes
**/

#define TP_VALID 0                /* density at touching point not zero */
#define TP_FZERO 1                /* denstiy 0 */
#define TP_GZERO 2                /* grandient of denstiy 0 */
#define TP_UGRAD 3                /* grandient of transformed denstiy does not fit into cone */
#define TP_NAN   9                /* numerical errors (NaN) */

#define TP_LEFT    1              /* minimum in left point */
#define TP_MIDDLE  2              /* minimum in middle point */
#define TP_RIGHT   3              /* minimum in right point */
#define TP_BRACKET 4              /* bracket found */
#define TP_EMPTY   0              /* no proper touching point found */
 
/*------------------------------------------------------------------*/

/** global variables                                               **/
/* /\* file handle for log file *\/ */
/* FILE *LOG; */

/* /\** debugging constanst                                            **\/ */
/* #define DB_VERTICES       2 */
/* #define DB_CONES          4 */
/* #define DB_EDGES          8 */
/* #define DB_CPARAMS        16 */
/* #define DB_RPOINT         32 */
/* #define DB_GUIDE          64 */
/* #define DB_GAMMA          128 */
/* #define DB_VOLUME         1024 */
/* #define DB_G              2048 */

/*------------------------------------------------------------------*/

/** Macros                                                         **/

/** Transformation                                                 **/
/*  this version supports T(x) = log(x) only                        */
#define T(x)       (log(x))       /* transformation function        */
#define T_deriv(x) (1./(x))       /* first derivative of transform funtion */
#define T_inv(x)   (exp(x))       /* inverse of transform function  */


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* since there is only file scope or program code, we abuse the              */
/* #include directive.                                                       */

/**  Public: User Interface (API)                                           **/
#include "mvtdr_newset.ch"

/**  Private                                                                **/
#include "mvtdr_init.ch"
#include "mvtdr_sample.ch"
#include "mvtdr_debug.ch"

/*---------------------------------------------------------------------------*/
