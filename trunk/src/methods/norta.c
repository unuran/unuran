/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      norta.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    generates random vector with given covariance matrix and     *
 *              given marginal distribution.                                 *
 *                                                                           *
 *   DESCRIPTION:                                                            *
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
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distributions/unur_distributions.h>
#include <uniform/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "norta.h"
/* #include "cstd.h" */
#include "auto.h"

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

#define NORTA_SET_NORMALGEN   0x001u    /* marginal normal variate generator */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NORTA"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_norta_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_norta_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_norta_nortu_setup( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Compute parameter for NORTU.                                              */
/*---------------------------------------------------------------------------*/

static void _unur_norta_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_norta_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_norta_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       par->data.norta       /* data for parameter object         */
#define GEN       gen->data.norta       /* data for generator object         */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

#define NORMAL    gen->gen_aux          /* pointer to normal variate generator */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_norta_new( const struct unur_distr *distr )
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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (!(distr->set & UNUR_DISTR_SET_RANKCORR)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"rank correlation matrix");
    return NULL; }

  /* TODO: */
  /*   if (!(distr->set & UNUR_DISTR_SET_MARGINAL)) { */
  /*     _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"marginals"); */
  /*     return NULL; } */

  /* allocate structure */
  par = _unur_xmalloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_NORTA_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR.normalgen = NULL;       /* use default marginal generator (=normal)    */

  par->method   = UNUR_METH_NORTA ;   /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_norta_init;

  return par;

} /* end of unur_norta_new() */

/*****************************************************************************/

struct unur_gen *
_unur_norta_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_NORTA ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NORTA_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_norta_create(par);
  if (!gen) { free(par); return NULL; }

  /* initialize normal generator if necessary */
  GEN.normaldistr = unur_distr_normal(NULL,0);
  if (NORMAL == NULL) {
    NORMAL = unur_init( unur_auto_new( GEN.normaldistr ) );
    if (NORMAL == NULL) {
      _unur_norta_free(gen);
      free(par);
      return NULL;
    }
  }
  /* need same uniform random number generator as NORTA generator */
  NORMAL->urng = gen->urng;
  /* copy debugging flags */
  NORMAL->debug = gen->debug;


  /* TODO: initialize generator for "marginal" distributions */

/*   /\* initialize generators for marginal distribution *\/ */
/*   if (_unur_distr_cvec_marginals_are_equal(DISTR.stdmarginals, GEN.dim)) { */
/*     /\* we can use the same generator object for all marginal distribuitons *\/ */
/*     /\** TODO: CSTD vor HINV vor NINV; nicht AUTO !! **\/  */
/*     struct unur_gen *marginalgen = unur_init( unur_auto_new( DISTR.stdmarginals[0] ) ); */
/*     if (marginalgen) */
/*       gen->gen_aux_list = _unur_gen_list_set(marginalgen,GEN.dim); */
/*   } */

/*   else { */
/*     int i,j; */
/*     int failed = FALSE; */
/*     struct unur_gen **marginalgens = _unur_xmalloc( GEN.dim * sizeof(struct unur_gen*) ); */
/*     for (i=0; i<GEN.dim; i++) { */
/*       /\** TODO: CSTD vor HINV vor NINV; nicht AUTO !! **\/  */
/*       marginalgens[i] = unur_init( unur_auto_new( DISTR.stdmarginals[i] ) ); */
/*       if (marginalgens[i]==NULL) { */
/*         failed=TRUE; break;  */
/*       } */
/*     } */
/*     if (failed) { */
/*       for (j=0; j<i; j++) _unur_free(marginalgens[j]); */
/*       free (marginalgens); */
/*     } */
/*     else */
/*       gen->gen_aux_list = marginalgens; */
/*   } */

/*   /\* the marginal generator is an auxiliary generator for method NORTA, of course *\/ */
/*   GEN.marginalgen_list = gen->gen_aux_list; */
  
/*   /\* verify initialization of marginal generators *\/ */
/*   if (GEN.marginalgen_list == NULL) { */
/*     _unur_error(gen->genid,UNUR_ERR_GENERIC,"init of marginal generators failed"); */
/*     _unur_norta_free(gen); */
/*     free(par); */
/*     return NULL; */
/*   } */

  /* compute parameters for NORTU (normal to uniform) */
  _unur_norta_nortu_setup( gen );


#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_norta_debug_init(gen);
#endif

  /* free parameters */
  free(par);

  /* o.k. */
  return gen;

} /* end of _unur_norta_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_norta_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NORTA_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_NORTA_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_norta_sample_cvec;
  gen->destroy = _unur_norta_free;
  gen->clone = _unur_norta_clone;

  /* dimension of distribution */
  GEN.dim = gen->distr->dim;

  /* copy data */
  NORMAL = PAR.normalgen;

  /* allocate array for auxiliary copula */
  GEN.copula = _unur_xmalloc(sizeof(double)*GEN.dim);

  /* initialize pointer */
  GEN.normaldistr = NULL;
  GEN.marginalgen_list = NULL;

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_norta_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_norta_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.norta

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NORTA_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate array for auxiliary copula */
  CLONE.copula = _unur_xmalloc(sizeof(double)*GEN.dim);

  /* clone marginal distribution */
  CLONE.normaldistr = _unur_distr_clone(GEN.normaldistr);



  /* cholesky factor of covariance matrix */
/*   CLONE.cholesky =  clone->distr->data.cvec.cholesky; */

  /* TODO: marginal gen and others !! */
  /* (normal generator) */

  /* marginal generators are (also) stored as auxiliary generator */
  /* which has already been cloned by generic_clone.              */
  CLONE.marginalgen_list = clone->gen_aux_list;

  return clone;

#undef CLONE
} /* end of _unur_norta_clone() */

/*---------------------------------------------------------------------------*/

int
_unur_norta_nortu_setup( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Compute parameter for NORTU (Normal to uniform)                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS if computed successfully                              */
     /*   error code   otherwise                                             */
     /*----------------------------------------------------------------------*/
{
  int dim = GEN.dim;    /* dimension of distribution */

  return UNUR_SUCCESS;
} /* end of _unur_norta_nortu_setup() */

/*****************************************************************************/

void
_unur_norta_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) (a*GEN.dim+b)
  int j;
  double *u;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  /* pointer to auxiliary array of uniforms */
  u = GEN.copula;

  /* generate random vector with independent components */
  for (j=0; j<GEN.dim; j++)
    u[j] = unur_sample_cont(NORMAL);


/*   /\* generate random vector with independent components *\/ */
/*   for (j=0; j<GEN.dim; j++) */
/*     vec[j] = unur_sample_cont(GEN.marginalgen_list[j]); */

/*   /\*  */
/*      transform to desired covariance structure:  */
/*      X = L.Y + mu  */
/*      where */
/*      L  ... cholesky factor of the covariance matrix */
/*      Y  ... vector with indenpent components (generated above) */
/*      mu ... mean vector */
/*      (notice that L is a lower triangular matrix) */
/*   *\/ */
/*   for (k=GEN.dim-1; k>=0; k--) { */
/*     vec[k] *= GEN.cholesky[idx(k,k)]; */
/*     for (j=k-1; j>=0; j--) */
/*       vec[k] += vec[j] * GEN.cholesky[idx(k,j)]; */
/*     vec[k] += DISTR.mean[k]; */
/*   } */

  for (j=0; j<GEN.dim; j++)
    vec[j] = u[j];

  return;

#undef idx
} /* end of _unur_norta_sample_cvec() */

/*****************************************************************************/

void
_unur_norta_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NORTA ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  /* free auxiliary arrays */
  if (GEN.copula) free (GEN.copula);

  /* free normal distribution object */
  if (GEN.normaldistr) free (GEN.normaldistr);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  _unur_generic_free(gen);

} /* end of _unur_norta_free() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_norta_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
/*   int i; */
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = NORTA (Vector Matrix Transformation)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

/*   fprintf(log,"%s: generators for standardized marginal distributions = \n",gen->genid); */
/*   fprintf(log,"%s:\t",gen->genid); */
/*   for (i=0; i<GEN.dim; i++) */
/*     fprintf(log,"[%s] ", GEN.marginalgen_list[i]->genid); */
/*   fprintf(log,"\n%s:\n",gen->genid); */

/*   fprintf(log,"%s: sampling routine = _unur_norta_sample()\n",gen->genid); */
/*   fprintf(log,"%s:\n",gen->genid); */

/*   fprintf(log,"%s:\n",gen->genid); */
/*   fprintf(log,"%s: INIT completed **********************\n",gen->genid); */

} /* end of _unur_norta_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/








#if 0

int
_unur_norta_nortu_setup( strunct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Compute parameter for NORTU.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS if computed successfully                              */
     /*   error code   otherwise                                             */
     /*----------------------------------------------------------------------*/
/* input 
 * dim ... dimension
 * cor ... correlation matrix (must be symmetric, positive definite and
 *         must have ones in the diagonal; for uniform marginals the 
 *         rank correlation and the (standard=Pearson) correlation are the same)
 * eps ... small positive constant, negative Eigenvalues are set to eps
 *
 * calculates the correlation matrix Sigma_Y for the normal vector such that
 * transformed to uniform the the correlation matrix cor is attained.
 * If Sigma_Y is not positive definite a "close" matrix is found by
 * calculating the Eigenvalues and Eigenvectors and forcing that all
 * Eigenvalues are non-negative. This idea is due to S. Ghosh but not published yet.
 *
 */
{ 
#define idx(a,b) ((a)*dim+(b))

  double *sigma_y;     /* covariance matrix for ... */

  int i,j,notposdef,notpos;
  double *helpm,*eval,*evec;
  struct nortu_gen *gen;
  
/*   int dim = GEN.dim; */

  /* allocate working arrays */
  sigma_y = _unur_xmalloc(dim*dim*sizeof(double));


/*   gen = malloc(sizeof(struct nortu_gen)); */
/*   gen->dim = dim; */
  gen->cholesky = malloc(sizeof(double)*dim*dim);
/*   gen->distrN01 = unur_distr_normal(NULL,0);/\* to use the CDF in nortu_sample*\/ */
/*   gen->genN01 = unur_str2gen("normal()");/\*to generate N(0,1) variates in nortu_sample*\/ */
  helpm = malloc(sizeof(double)*dim*dim);

  //  printf("begin of nortusetup correlation matrix:\n");
  //  printsqmatrix(cor,dim);

  /* setup correlation matrix for ... */
  for(i=0; i<dim; i++) {
    for(j=0; j<i; j++)
      /* left lower part: make matrix symmetric */
	 sigma_y[idx(i,j)] = sigma_y[idx(j,i)];
    /*   diagonal */
    sigma_y[idx(i,i)] = 1.;
    for(j=i+1;j<dim;j++)
      /* right upper part */
      sigma_y[idx(i,j)] = 2.*sin(GEN.corr[idx(i,j)]*(M_PI/6.));  
  }


  //  printf("nortusetup sigma_y:\n");
  //  printsqmatrix(sigma_y,dim);

  gsl_matrix_view m
    = gsl_matrix_view_array (sigma_y, dim, dim);
  gsl_vector *evalgsl = gsl_vector_alloc (dim);
  gsl_matrix *evecgsl = gsl_matrix_alloc (dim, dim);
  gsl_eigen_symmv_workspace * w =
    gsl_eigen_symmv_alloc (dim);
  gsl_eigen_symmv (&m.matrix, evalgsl, evecgsl, w);
  gsl_eigen_symmv_free (w);
  eval = evalgsl->data;
  evec = evecgsl->data;

  //  printf("the eigenvalues\n");
  //  printvector(eval,dim);

  /*  printf("the eigenvectors\n");
  printsqmatrix(evec,dim);
  */
  /* the vector eval is now containing the eigenvalues
     and evec the corresponding eigenvectors */

  /*as gsl_eigen_symmv is destroying the diagonal and the lower triangular
    we have to restore the full matrix sigma_y  */
  for(i=0;i<dim;i++){
    for(j=0;j<i;j++)
      sigma_y[idx(i,j)] = sigma_y[idx(j,i)];
    sigma_y[idx(i,i)] = 1.;
  }


  /* check if all eval are positive */
  notpos = 0;
  for(i=0; i<dim; i++)
    if(eval[i]<=0.){ /*oder <= eps ?????*/
      eval[i]=eps;
      notpos++;
    }

  if(notpos>0){
    mmult_matr_diagonal(evec,eval,dim,helpm);
    mmult_matr_matrt(helpm,evec,dim,sigma_y);
    makecorrelationmatrix(sigma_y,dim);
    //    printf("as sigma_y not positive definite EC corrected correlation matrix\n");
    //printsqmatrix(sigma_y,dim);
  }


  //  printf("nortusetup sigma_y:\n");
  //  printsqmatrix(sigma_y,dim);
  notposdef = cholesky_decomp(sigma_y, dim, gen->cholesky);
  //  printf("nortusetup cholesky factor:\n");
  //  printsqmatrix(gen->cholesky,dim);



  gsl_vector_free(evalgsl);
  gsl_matrix_free(evecgsl);
  free(sigma_y);
  free(helpm);
  if(notposdef){
    printf("ERROR nortu-setup: Eigenvalue corrected matrix not positive definite\n");
    free(gen->cholesky);
    free_nortu_gen(gen);
    return NULL;
  }
  return gen;
#undef idx
} /* end of _unur_norta_nortu_setup() */

/*****************************************************************************/

void
_unur_norta_sample_cvec( struct unur_gen *gen, double *vec )
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
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  /* generate random vector with independent components */
  for (j=0; j<GEN.dim; j++)
    vec[j] = unur_sample_cont(GEN.uvgen);

  /* 
     transform to desired covariance structure: 
     X = L.Y + mu 
     where
     L  ... cholesky factor of the covariance matrix
     Y  ... vector with indenpent components (generated above)
     mu ... mean vector
     (notice that L is a lower triangular matrix)
  */
  for (k=GEN.dim-1; k>=0; k--) {
    vec[k] *= GEN.cholesky[idx(k,k)];
    for (j=k-1; j>=0; j--)
      vec[k] += vec[j] * GEN.cholesky[idx(k,j)];
    vec[k] += DISTR.mean[k];
  }

#undef idx
} /* end of _unur_norta_sample_cvec() */

/*****************************************************************************/

#endif
