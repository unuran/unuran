/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vmt.h                                                        *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    generated random vecto with independent components with      *
 *              given marginal distribution and use linear transformation    *
 *              of vector.                                                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      multivariate distribution with given mean vector and                 *
 *      covariance matrix.                                                   *
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
 *   [1] Devroye, L. (1986): Non-Uniform Random Variate Generation, New-York *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

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

#define VMT_SET_UVGEN        0x001u    /* marginal generator set             */

/*---------------------------------------------------------------------------*/

#define GENTYPE "VMT"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vmt_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vmt_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vmt_default_uvgen( void );
/*---------------------------------------------------------------------------*/
/* get default generator for "marginal" distribution.                        */
/*---------------------------------------------------------------------------*/

static void _unur_vmt_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_vmt_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double *cholesky_decomposition( double *S, int dim );
/*---------------------------------------------------------------------------*/
/* the Colesky factor of a covariance matrix S is computed and returned      */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_vmt_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       par->data.vmt         /* data for parameter object         */
#define GEN       gen->data.vmt         /* data for generator object         */
#define DISTR     gen->distr.data.cvec  /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_vmt_new( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
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
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (!(distr->set & UNUR_DISTR_SET_MEAN)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mean"); return NULL; }
  if (!(distr->set & UNUR_DISTR_SET_COVAR)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"covariance matrix");
    return NULL; }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_VMT_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR.uvgen     = NULL;       /* use default marginal generator (=normal)    */

  par->method   = UNUR_METH_VMT ;     /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_vmt_init;

  return par;

} /* end of unur_vmt_new() */

/*****************************************************************************/

int
unur_vmt_set_marginalgen( struct unur_par *par, struct unur_gen *uvgen )
     /*----------------------------------------------------------------------*/
     /* set generator for (univariate) marginal distribution.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   uvgen ... pointer to generator for (univariate) marginal distr.    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,VMT );

  /* check new parameter for generator */

  if ( (uvgen->method & UNUR_MASK_TYPE) != UNUR_METH_CONT ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"marginal generator");
    return 0;
  }
    
  /* store date */
  PAR.uvgen = uvgen;

  /* changelog */
  par->set |= VMT_SET_UVGEN;

  return 1;

} /* end of unur_vmt_set_marginalgen() */

/*****************************************************************************/

struct unur_gen *
_unur_vmt_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_VMT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VMT_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_vmt_create(par);
  if (!gen) { free(par); return NULL; }

  /* initialize generator for "marginal" distribution */
  if (GEN.uvgen == NULL) {
    GEN.uvgen = _unur_vmt_default_uvgen();
    if (GEN.uvgen == NULL) {
      _unur_error(gen->genid,UNUR_ERR_GENERIC,"init of marginal generator failt");
      unur_free(gen);
      return NULL;
    }
  }
  /* else: generator provided by user */

  /* cholesky factor of covariance matrix */
  if (DISTR.covar)
    GEN.cholesky = cholesky_decomposition( DISTR.covar, GEN.dim );
      
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_vmt_debug_init(gen);
#endif

  /* free parameters */
  free(par);

  if (DISTR.covar && !GEN.cholesky) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"covariance matrix not positive definite");
    _unur_vmt_free(gen);
    return NULL;
  }

  /* o.k. */
  return gen;

} /* end of _unur_vmt_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_vmt_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VMT_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_VMT_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* dimension of distribution */
  GEN.dim = gen->distr.dim; 

  /* mean vector */
  if (DISTR.mean) {
    DISTR.mean = _unur_malloc( GEN.dim * sizeof(double) );
    memcpy( DISTR.mean, par->distr->data.cvec.mean, GEN.dim * sizeof(double) );
  }
  /* else: DISTR.mean == NULL  --> (0,0,...,0) */

  if (DISTR.covar) {
    DISTR.covar = _unur_malloc( GEN.dim * GEN.dim * sizeof(double) );
    memcpy( DISTR.covar, par->distr->data.cvec.covar, GEN.dim * GEN.dim * sizeof(double) );
  }
  /* else: DISTR.covar == NULL  --> identity matrix */

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_vmt_sample_cvec;
  gen->destroy = _unur_vmt_free;

  /* copy some parameters into generator object */
  GEN.uvgen = PAR.uvgen;            /* generator for univariate distribution */

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* mean vector for generator */
  GEN.mean = DISTR.mean;

  /* initialize pointer */
  GEN.cholesky = NULL;

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_vmt_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vmt_default_uvgen( void )
     /*----------------------------------------------------------------------*/
     /* initialize generator object for "marginal" distribution              */
     /* default is standard normal distribution                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   none                                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object (standard normal)                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_DISTR *distr;
  UNUR_GEN *uvgen;

  /* construct a generator object for univariate distribuiton */
  distr = unur_distr_normal(NULL,0);
  uvgen = unur_init( unur_cstd_new( distr ) );
  unur_distr_free(distr);

  /* o.k. */
  return uvgen;
} /* end of _unur_vmt_default_uvgen() */

/*****************************************************************************/

void
_unur_vmt_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) (a*GEN.dim+b)
  int j,k;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_VMT_GEN,/*void*/);

  /* generate random vector with independent components */
  for (j=0; j<GEN.dim; j++)
    vec[j] = unur_sample_cont(GEN.uvgen);

  /* transform to desired covariance structure */
  /* notice that GEN.cholesky is a lower triangular matrix */
  for (k=GEN.dim-1; k>=0; k--) {
    vec[k] *= GEN.cholesky[idx(k,k)];
    for (j=k-1; j>=0; j--)
      vec[k] += vec[j] * GEN.cholesky[idx(k,j)];
    vec[k] += GEN.mean[k];
  }

#undef idx
} /* end of _unur_vmt_sample_cvec() */

/*****************************************************************************/

void
_unur_vmt_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_VMT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VMT_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (!(gen->set & VMT_SET_UVGEN)) unur_free(GEN.uvgen);
  if (DISTR.mean) free(DISTR.mean);
  if (DISTR.covar) free(DISTR.covar);
  free(GEN.cholesky);
  _unur_free_genid(gen);
  free(gen);

} /* end of _unur_vmt_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

double *
cholesky_decomposition( double *S, int dim )
     /*----------------------------------------------------------------------*/
     /* the Colesky factor L of a variance-covariance matrix S is computed   */
     /* (the necessary array is allocated)                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   S   ... variance-covariance matrix                                 */
     /*   dim ... dimension (S is a dim x dim matrixes)                      */ 
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to chokesky factor                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) (a*dim+b)

  double *L;
  int i,j,k;
  double sum1,sum2;

  /* allocate memory for cholesky factor */
  L = _unur_malloc( dim * dim * sizeof(double) );

  /* run cholesky decomposition */
  L[idx(0,0)] = sqrt( S[idx(0,0)] );

  for(j=1; j<dim; j++) {
    L[idx(j,0)] = S[idx(j,0)] / L[idx(0,0)];

    sum1 = L[idx(j,0)] * L[idx(j,0)];
    for(k=1; k<j; k++) {
      sum2 = 0.;
      for(i=0; i<k; i++)
	sum2 += L[idx(j,i)] * L[idx(k,i)];
      
      L[idx(j,k)] = (S[idx(j,k)] - sum2) / L[idx(k,k)];
      sum1 += L[idx(j,k)] * L[idx(j,k)];
    }

    if (S[idx(j,j)] <= sum1) {
      /* covariance matrix not positive definite */
      free(L); return NULL;
    }

    L[idx(j,j)] = sqrt( S[idx(j,j)] - sum1 );
  }

  /* although not necessary upper triangular of L - matrix is set to 0 */
  for(j=0; j<dim; j++)
    for(k=j+1; k<dim; k++)
      L[idx(j,k)]=0.;

  /* return (pointer to) cholesky factor */
  return L;

#undef idx
} /* end of cholesky_decomposition() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_vmt_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_VMT_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = VMT (Vector Matrix Transformation)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s:\tcholesky factor = ",gen->genid);
  if (DISTR.covar && !GEN.cholesky) {
    fprintf(log,"[matrix not positive definite !]\n");
  }
  else {
    fprintf(log,"%s\n",(GEN.cholesky ? "" : "[NULL]"));
    for (j=0; j<GEN.dim; j++) {
      fprintf(log,"%s:\t   (%7.4f",gen->genid,(GEN.cholesky ? GEN.cholesky[GEN.dim*j] : (j==0?1.:0.)));
      for (i=1; i<GEN.dim; i++) 
	fprintf(log,",%7.4f",(GEN.cholesky ? GEN.cholesky[GEN.dim*j+i] : (i==j?1.:0.)));
      fprintf(log,")\n");
    }
  }
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: marginal distribution = %s    [genid = %s]\n",gen->genid,GEN.uvgen->distr.name,GEN.uvgen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = _unur_vmt_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: INIT completed **********************\n",gen->genid);

} /* end of _unur_vmt_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
