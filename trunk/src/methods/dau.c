/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dau.h                                                        *
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
 *   PARAMETERS:                                                             *
 *      int factor  ... relative size of alias table (relative to N)         *
 *                      (default: 1)                                         *
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
 *                                                
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Wed Sep  8 16:15:41 CEST 1999                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

#define GENTYPE "DAU"         /* type of generator                           */

/*---------------------------------------------------------------------------*/

static struct unur_gen *dau_create( struct unur_par *par);
#if UNUR_DEBUG & UNUR_DB_INFO
static void dau_info_init( struct unur_par *par, struct unur_gen *gen );
static void dau_info_table( struct unur_gen *gen );
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR     par->data.dau
#define GEN     gen->data.dau
#define SAMPLE  gen->sample.discr

/*---------------------------------------------------------------------------*/
/* Special debugging flags (do not use the first 3 bits) */
#define DAU_DB_TABLE   0x010UL

#define HIST_WIDTH   40  /* width of histogram for printing alias table      */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
 unur_dau_new( double *prob, int len )
/*---------------------------------------------------------------------------*/
/* get default parameters                                                    */
/*                                                                           */
/* parameters:                                                               */
/*   prob ... pointer to probability vector                                  */
/*   len  ... length of probability vector                                   */
/*                                                                           */
/* return:                                                                   */
/*   default parameters (pointer to structure)                               */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_DAU_PAR);

  /* copy input */
  PAR.prob       = prob;
  PAR.len        = len;

  /* set default values */
  PAR.urn_factor = 1.;              /* use same size for table               */

  par->method    = UNUR_METH_DAU;   /* method and default variant            */
  par->set       = 0UL;             /* inidicate default parameters          */    
  par->urng      = unur_get_default_urng(); /* use default urng              */

  _unur_set_debugflag_default(par); /* set default debugging flags           */
  _unur_set_genid(par,GENTYPE);     /* set generator identifier              */

  /* routine for starting generator */
  par->init = unur_dau_init;

  return par;

} /* end of unur_dau_new() */

/*****************************************************************************/

struct unur_gen *
 unur_dau_init( struct unur_par *par )
/*---------------------------------------------------------------------------*/
/* initialize new generator                                                  */
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*                                                                           */
/* return:                                                                   */
/*   pointer to generator object                                             */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  int *begin, *poor, *rich;   /* list of (rich and poor) strips */
  int *npoor;                   /* next poor on stack */
  double sum, ratio;
  int i;                      /* aux variable */

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_DAU_PAR,NULL);

  /* create a new empty generator object */
  gen = dau_create(par);
  if (!gen) { free(par); return NULL; }

  /* compute sum of all probabilities */
  for( sum=0, i=0; i<GEN.len; i++ ) {
    sum += PAR.prob[i];
#if CHECKARGS
    /* ... and check probability vector */
    if (PAR.prob[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"probability < 0 not possible!");
      free(par); unur_dau_free(gen);
      return NULL;
    }
#endif
  }

  /* make list of poor and rich strips */
  begin = _unur_malloc( (GEN.urn_size+2) * sizeof(int) );
  poor = begin;                    /* poor strips are stored at the beginning ... */
  rich = begin + GEN.urn_size + 1; /* ... rich strips at the end of the list      */

  /* copy probability vector; scale so that it sums to GEN.urn_size and */
  /* find rich and poor strips at start                                 */
  ratio = GEN.urn_size / sum;
  for( i=0; i<GEN.len; i++ ) {
    GEN.qx[i] = PAR.prob[i] * ratio;
    if (GEN.qx[i] >= 1.) {  /* rich strip                  */
      *rich = i;            /* add to list ...             */
      --rich;               /* and update pointer          */
      GEN.jx[i] = i;        /* init donor (itself)         */
    }
    else {                  /* poor strip                  */
      *poor = i;            /* add to list                 */
      ++poor;               /* update pointer              */
      /* it is not necessary to mark donor                 */
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
    /* this must not happen */
    _unur_error(gen->genid,UNUR_ERR_INIT,"no rich strips found for Robin Hood algorithm.");
    unur_dau_free(gen); free(par); free(begin);
    return NULL;
  }
  
  /* rich doesn not point to a valid rich strip yet */
  ++rich;        /* now it does                     */

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
    if (fabs(sum) > TOLERANCE) {
      /* sum of deviations too large --> serious error */
      _unur_error(gen->genid,UNUR_ERR_INIT,"serious roundoff error when making squared histogram!");
      unur_dau_free(gen); free(par); free(begin);
      return NULL; 
    }
  }

  /* free lists of strips */
  free(begin);
  
  /* write info into log file */
#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) dau_info_init(par,gen);
#endif

  /* free parameters */
  free(par);
  
  return gen;

} /* end of unur_dau_init() */

/*****************************************************************************/

int
 unur_dau_sample( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*                                                                           */
/* return:                                                                   */
/*   integer (sample from random variate)                                    */
/*                                                                           */
/* error:                                                                    */
/*   return 0                                                                */
/*---------------------------------------------------------------------------*/
{ 
  int iu;
  double u;

  /* check arguments */
  CHECK_NULL(gen,0);
  COOKIE_CHECK(gen,CK_DAU_GEN,0);

  /* sample from U(0,urn_size) */
  u = _unur_call_urng(gen);
  u *= GEN.urn_size;
  iu = (int) u;

  /* immediate return ? */
  if (iu >= GEN.len) return (GEN.jx[iu]);

  /* else choose number or its alias at random */
  u -= iu;   /* reuse of random number */
  return ((u <= GEN.qx[iu]) ? iu : GEN.jx[iu] );

} /* end of unur_dau_sample() */

/*****************************************************************************/

void
 unur_dau_free( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* deallocate generator object                                               */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*---------------------------------------------------------------------------*/
{ 

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* magic cookies */
  COOKIE_CHECK(gen,CK_DAU_GEN,/*void*/);
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(GEN.jx);
  free(GEN.qx);
  if (GEN.prob) free(GEN.prob);
  free(gen);

} /* end of unur_dau_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
 dau_create( struct unur_par *par)
/*---------------------------------------------------------------------------*/
/* allocate memory for generator                                             */
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*                                                                           */
/* return:                                                                   */
/*   pointer to (empty) generator object with default settings               */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_DAU_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc(sizeof(struct unur_gen));

  /* magic cookies */
  COOKIE_SET(gen,CK_DAU_GEN);

  /* copy some parameters into generator object */
  GEN.len = PAR.len;                /* length of probability vector               */
  GEN.prob = NULL;                  /* copy probability vector on demand          */
  gen->method = par->method;        /* indicates used method                      */
  _unur_copy_urng_pointer(par,gen); /* copy pointer to urng into generator object */
  _unur_copy_debugflag(par,gen);    /* copy debugging flags into generator object */
  _unur_copy_genid(par,gen);        /* copy generator identifier                  */

  /* routines for sampling and destroying generator */
  SAMPLE = unur_dau_sample;
  gen->destroy = unur_dau_free;

  /* size of table */
  GEN.urn_size = (int)(PAR.len * PAR.urn_factor);
  if (GEN.urn_size < PAR.len)
    /* do not use a table that is smaller then length of probability vector */
    GEN.urn_size = PAR.len;

  /* allocate memory for the tables */
  GEN.jx = _unur_malloc( GEN.urn_size * sizeof(int) );
  GEN.qx = _unur_malloc( GEN.urn_size * sizeof(double) );

  /* copy probability vector on demand */
  if (par->method & UNUR_MASK_COPYALL) {
    /* need an arry for the probabilty vector */
    GEN.prob = _unur_malloc( GEN.len * sizeof(double) );
    /* copy */
    for (i=0; i<GEN.len; i++)
      GEN.prob[i] = PAR.prob[i];
  }
  
  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of dau_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

static void
 dau_info_init( struct unur_par *par, struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* write info about generator into logfile                                   */
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*   gen ... pointer to generator object                                     */
/*---------------------------------------------------------------------------*/
{
  FILE *log;

  log = unur_get_log();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variate\n",gen->genid);
  fprintf(log,"%s: method  = alias and alias-urn method\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = unur_dau_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: length of probability vector = %d\n",gen->genid,PAR.len);
  fprintf(log,"%s: size of urn table = %d   (rel. = %g%%",
	  gen->genid,GEN.urn_size,100.*PAR.urn_factor);
  _unur_print_if_default(par,UNUR_SET_FACTOR);
  if (GEN.urn_size == GEN.len)
    fprintf(log,")   (--> alias method)\n");
  else
    fprintf(log,")   (--> alias-urn method)\n");
  fprintf(log,"%s:\n",gen->genid);

  if (gen->debug & DAU_DB_TABLE) {
    dau_info_table(gen);
    fprintf(log,"%s:\n",gen->genid);
  }

} /* end of dau_info_init() */

/*---------------------------------------------------------------------------*/

static void
 dau_info_table( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* print alias table into logfile                                           */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*---------------------------------------------------------------------------*/
{
  FILE *log;
   
  int i, j, m;
  
  log = unur_get_log();
  
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

} /* end of dau_info_table() */

#endif

/*****************************************************************************/
