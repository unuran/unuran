/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dau.c                                                        *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    alias and alias-urn method                                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given N discrete events with different probabilities P[k]            *
 *      produce a value k consistent with its probability.                   *
 *                                                                           *
 *   REQUIRED:  pointer to probability vector                                *
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
 *   [1] Walker, A. J. (1974): New fast method for generating discrete       *
 *       random numbers with arbitrary frequency distributions,              *
 *       Electron. Lett. 10, pp. 127-128                                     *
 *                                                                           *
 *   [2] Walker, A. J. (1977): An efficient method for generating discrete   *
 *       random variables with general distributions,                        *
 *       ACM Trans. Math. Software 3, pp. 253-256                            *
 *                                                                           *
 *   [3] Kronmal, R. A. and Peterson, A. V. (1979): On the alias method for  *
 *       generating random variables from a discrete distribution,           *
 *       Amer. Statist. 33(4), pp. 214-218                                   *
 *                                                                           *
 *   [4] Peterson, A. V. and Kronmal, R. A. (1982): On mixture methods for   *
 *       the computer generation of random variables,                        *
 *       Amer. Statist. 36(3), pp. 184-191                                   *
 *                                                                           *
 *   SUMMARY:                                                                *
 *   [5] Devroye, L. (1986): Non-Uniform Random Variate Generation, New-York *
 *                                                                           *
 *   [6] Knuth, D. E. (1997): The Art of Computer Programming,               *
 *       Volume 2 (Seminumerical algorithms), 3rd edition,                   *
 *       Addison-Wesley (1997), p.120                                        *
 *                                                                           *
 *   SEE ALSO:                                                               *
 *   [7] Zaman, A. (1996), Generation of Random Numbers from an Arbitrary    *
 *       Unimodal Density by Cutting Corners, unpublished manuskript         *
 *       available at http://chenab.lums.edu.pk/~arifz/                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * This algorithmus is based on [1,2]. It is an ingeneous method for         *
 * generating random variates with finite probability vector which requires  *
 * a table of size N and needs only one comparision.                         *
 *                                                                           *
 * Walker's algorithm does some preprocessing, and provides two array:       *
 * floating point Q[k] and integer J[k]. A value k is chosen from 0..N-1     *
 * with equal likelihood, and then a uniform random number u is compared to  *
 * Q[k].  If it is less than Q[k], then k is returned. Otherwise, J[k] is    *
 * returned.                                                                 *
 * The method has been generalized in [4] to the "alias-urn" method that     *
 * works in the exactly same way but uses a larger table of aliases.         *
 *                                                                           *
 * The original algorithm needs 2 uniform random numbers. By reusing only    *
 * one is necessary (see [6]).                                               *
 *                                                                           *
 * Walker's original paper describes an O(N^2) algorithm for setting         *
 * up the Q and J arrays. It had been improved in e.g. [3].                  *
 * This implementation uses Marsaglia's "Robin Hood algorithm" (see [7]):    *
 * This O(N) algorithm goes through all the p_k's and decides if they are    *
 * are "poor" or "rich" according to whether they are less than or greater   *
 * than (>=) the mean value 1 (For convienience we normalize the p_k's,      *
 * s.t. there sum N instead of 1.). The indices to the poors and the richs   *
 * are put in separate stacks, and then we work through the stacks together. *
 * We take from the next rich on the stack and give to the poor, s.t. it     *
 * has the average value, and then it is popped from the stack of poors and  *
 * stored in the tables Q and J.                                             *
 * (Walker always wanted to pair up the poorest with the richest.)           *
 * This reduces the size of the rich and even might become poor, i.e.,       *
 * it is popped from the stack of richs and pushed on the stack of poors.    *
 * Since the size of the the two stacks together is <= N, we use one array   *
 * and store the poors on the beginning and the richs at the and of an       *
 * array of size N.                                                          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <uniform/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dau.h"

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DAU_DEBUG_PRINTVECTOR   0x00000100u
#define DAU_DEBUG_TABLE         0x00000200u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DAU_SET_URNFACTOR       0x01u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DAU"         /* type of generator                           */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dau_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dau_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dau_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_dau_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dau_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_dau_debug_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print data for alias table.                                               */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr      /* data for distribution object      */

#define PAR       par->data.dau         /* data for parameter object         */
#define GEN       gen->data.dau         /* data for generator object         */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dau_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL(GENTYPE,distr,NULL);

  /* check distribution */
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);

  if (DISTR_IN.pv == NULL) {
    /* There is no PV try to compute it.                         */
    if ( DISTR_IN.pmf
	 && ( (((unsigned)DISTR_IN.domain[1] - (unsigned)DISTR_IN.domain[0]) < UNUR_MAX_AUTO_PV)
	      || ( (distr->set & UNUR_DISTR_SET_PMFSUM) && DISTR_IN.domain[0] > INT_MIN ) ) ) {
      /* However this requires a PMF and either a bounded domain   */
      /* or the sum over the PMF.                                  */
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV. Try to compute it.");
    }
    else {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); return NULL;
    }
  }

  /* allocate structure */
  par = _unur_xmalloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_DAU_PAR);

  /* copy input */
  par->distr     = distr;            /* pointer to distribution object       */

  /* set default values */
  PAR.urn_factor = 1.;               /* use same size for table              */

  par->method    = UNUR_METH_DAU;    /* method                               */
  par->variant   = 0u;               /* default variant (no other variants)  */
  par->set       = 0u;               /* inidicate default parameters         */    
  par->urng      = unur_get_default_urng(); /* use default urng              */
  par->urng_aux  = NULL;                    /* no auxilliary URNG required   */

  par->debug     = _unur_default_debugflag; /* set default debugging flags   */

  /* routine for starting generator */
  par->init = _unur_dau_init;

  return par;

} /* end of unur_dau_new() */

/*---------------------------------------------------------------------------*/

int
unur_dau_set_urnfactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of urn                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of urn                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DAU );
  
  /* check new parameter for generator */
  if (factor < 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"relative urn size < 1.");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR.urn_factor = factor;

  /* changelog */
  par->set |= DAU_SET_URNFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_dau_set_urnfactor() */

/*****************************************************************************/

struct unur_gen *
_unur_dau_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  int *begin, *poor, *rich;     /* list of (rich and poor) strips */
  int *npoor;                   /* next poor on stack */
  double *pv;                   /* pointer to probability vector */
  int n_pv;                     /* length of probability vector */
  double sum, ratio;
  int i;                        /* aux variable */

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_DAU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DAU_PAR,NULL);
  
  /* create a new empty generator object */
  gen = _unur_dau_create(par);
  if (!gen) { free(par); return NULL; }

  /* probability vector */
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;

  /* compute sum of all probabilities */
  for( sum=0, i=0; i<n_pv; i++ ) {
    sum += pv[i];
    /* ... and check probability vector */
    if (pv[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"probability < 0");
      free(par); _unur_dau_free(gen);
      return NULL;
    }
  }

  /* make list of poor and rich strips */
  begin = _unur_xmalloc( (GEN.urn_size+2) * sizeof(int) );
  poor = begin;                    /* poor strips are stored at the beginning ... */
  rich = begin + GEN.urn_size + 1; /* ... rich strips at the end of the list      */

  /* copy probability vector; scale so that it sums to GEN.urn_size and */
  /* find rich and poor strips at start                                 */
  ratio = GEN.urn_size / sum;
  for( i=0; i<n_pv; i++ ) {
    GEN.qx[i] = pv[i] * ratio;  /* probability rescaled        */
    if (GEN.qx[i] >= 1.) {      /* rich strip                  */
      *rich = i;                /* add to list ...             */
      --rich;                   /* and update pointer          */
      GEN.jx[i] = i;            /* init donor (itself)           */
    }
    else {                      /* poor strip                    */
      *poor = i;                /* add to list                 */
      ++poor;                   /* update pointer              */
      /* it is not necessary to mark donor                     */
    }
  }

  /* all other (additional) strips own nothing yet */
  for( ; i<GEN.urn_size; i++ ) {
    GEN.qx[i] = 0.;
    *poor = i; 
    ++poor;
  }

  /* there must be at least one rich strip */
  if (rich == begin + GEN.urn_size + 1 ) {
    /* this must not happen:
       no rich strips found for Robin Hood algorithm. */
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_dau_free(gen); free(par); free(begin);
    return NULL;
  }
  
  /* rich must point to the first rich strip yet */
  ++rich;

  /* now make the "squared histogram" with Robin Hood algorithm (Marsaglia) */
  while (poor != begin) {
    if (rich > begin + GEN.urn_size + 1) {
      /* there might be something wrong; assume a neglectable round off error */
      break;
    }

    npoor = poor - 1;                       /* take next poor from stack */
    GEN.jx[*npoor] = *rich;                 /* store the donor */
    GEN.qx[*rich] -= 1. - GEN.qx[*npoor];   /* update rich */

    /* rich might has given too much, so it is poor then */
    if (GEN.qx[*rich] < 1.) {
      *npoor = *rich;      /* exchange noveau-poor with former poor in list */
      ++rich;              /* remove it from list of rich */
    }
    else
      --poor;              /* remove poor from list */
  }

  /* if there has been an round off error, we have to complete the table */
  if (poor != begin) {
    sum = 0.;                   /* we estimate the round off error            */
    while (poor != begin) {
      npoor = poor - 1;         /* take next poor from stack */
      sum += 1. - GEN.qx[*npoor];
      GEN.jx[*npoor] = *npoor;  /* mark donor as "not valid" */
      GEN.qx[*npoor] = 1.;      /* set probability to 1 (we assume that it is very close to one) */
      --poor;                   /* remove from list */
    }
    if (fabs(sum) > UNUR_SQRT_DBL_EPSILON)
      /* sum of deviations very large */
      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"squared histogram");
  }

  /* free lists of strips */
  free(begin);
  
  /* write info into log file */
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_dau_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);
  
  return gen;

} /* end of _unur_dau_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_dau_create( struct unur_par *par)
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DAU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_DAU_GEN);

  /* we need a PV */
  if (DISTR.pv == NULL) {
    /* try to compute PV */
    if (unur_distr_discr_make_pv( gen->distr ) <= 0) {
      /* not successful */
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); 
      _unur_distr_free(gen->distr); free(gen);
      return NULL;
    }
  }

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_dau_sample;
  gen->destroy = _unur_dau_free;
  gen->clone = _unur_dau_clone;

  /* copy some parameters into generator object */
  GEN.len = DISTR.n_pv;             /* length of probability vector          */

  /* size of table */
  GEN.urn_size = (int)(GEN.len * PAR.urn_factor);
  if (GEN.urn_size < GEN.len)
    /* do not use a table that is smaller then length of probability vector */
    GEN.urn_size = GEN.len;

  /* allocate memory for the tables */
  GEN.jx = _unur_xmalloc( GEN.urn_size * sizeof(int) );
  GEN.qx = _unur_xmalloc( GEN.urn_size * sizeof(double) );

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_dau_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dau_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE clone->data.dau

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DAU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy data for generator */
  CLONE.jx = _unur_xmalloc( GEN.urn_size * sizeof(int) );
  memcpy( CLONE.jx, GEN.jx, GEN.urn_size * sizeof(int) );
  CLONE.qx = _unur_xmalloc( GEN.urn_size * sizeof(double) );
  memcpy( CLONE.qx, GEN.qx, GEN.urn_size * sizeof(double) );

  return clone;

#undef CLONE
} /* end of _unur_dau_clone() */

/*****************************************************************************/

int
_unur_dau_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   integer (sample from random variate)                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{ 
  int iu;
  double u;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DAU_GEN,INT_MAX);

  /* sample from U(0,urn_size) */
  u = _unur_call_urng(gen->urng);
  u *= GEN.urn_size;
  iu = (int) u;

  /* immediate return ? */
  if (iu >= GEN.len) return (GEN.jx[iu] + DISTR.domain[0]);

  /* else choose number or its alias at random */
  u -= iu;   /* reuse of random number */

  return (((u <= GEN.qx[iu]) ? iu : GEN.jx[iu] ) + DISTR.domain[0]);

} /* end of _unur_dau_sample() */

/*****************************************************************************/

void
_unur_dau_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_DAU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DAU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free two auxiliary tables */
  if (GEN.jx) free(GEN.jx);
  if (GEN.qx) free(GEN.qx);

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_dau_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_dau_debug_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_DAU_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DAU_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variate\n",gen->genid);
  fprintf(log,"%s: method  = alias and alias-urn method\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_discr_debug( gen->distr,gen->genid,(gen->debug & DAU_DEBUG_PRINTVECTOR));

  fprintf(log,"%s: sampling routine = _unur_dau_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: length of probability vector = %d\n",gen->genid,GEN.len);
  fprintf(log,"%s: size of urn table = %d   (rel. = %g%%",
	  gen->genid,GEN.urn_size,100.*PAR.urn_factor);
  _unur_print_if_default(par,DAU_SET_URNFACTOR);
  if (GEN.urn_size == GEN.len)
    fprintf(log,")   (--> alias method)\n");
  else
    fprintf(log,")   (--> alias-urn method)\n");
  fprintf(log,"%s:\n",gen->genid);

  if (gen->debug & DAU_DEBUG_TABLE) {
    _unur_dau_debug_table(gen);
    fprintf(log,"%s:\n",gen->genid);
  }

} /* end of _unur_dau_debug_init() */

/*---------------------------------------------------------------------------*/

#define HIST_WIDTH   40  /* width of histogram for printing alias table      */

static void
_unur_dau_debug_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print alias table into logfile                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i, j, m;
  
  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DAU_GEN,RETURN_VOID);
   
  log = unur_get_stream();
  
  /* print head of table */
  fprintf(log,"%s: alias table:\n", gen->genid); 
  fprintf(log,"%s:\n", gen->genid);
  fprintf(log,"%s:         ratio donor/acceptor",gen->genid);
  for (i=0; i<HIST_WIDTH-17; i++)
    /* print blanks */
    fprintf(log," ");
  fprintf(log,"jx:     qx:\n");
    
  for (i=0; i<GEN.urn_size; i++){
    m = HIST_WIDTH * GEN.qx[i] + 0.5;
    fprintf(log,"%s:[%4d]: ", gen->genid,i); 

    /* illustrate ratio donor/acceptor graphically */
    for (j=0; j<HIST_WIDTH; j++)
      if (j<m)
	fprintf(log, "*"); 
      else                
	fprintf(log,"-");
 
    fprintf(log," %5d  ", GEN.jx[i]);           /* name donor */
    fprintf(log,"  %6.3f%%\n", GEN.qx[i]*100);  /* cut point */
  }

} /* end of _unur_dau_debug_table() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
