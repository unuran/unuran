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
#include <distr/cont.h>
#include <distributions/unur_distributions.h>
#include <uniform/urng.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "norta.h"
#include "vmt.h"

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* smallest eigenvalue allowed for correlation matrix                        */
#define UNUR_NORTA_MIN_EIGENVALUE  (1.e-10)

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define NORTA_DEBUG_SIGMA_Y     0x00000010u   /* print sigma_y for normal    */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

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
/* Compute parameter for NORTU (normal to uniform).                          */
/*---------------------------------------------------------------------------*/

static int _unur_norta_make_correlationmatrix( int dim, double *M);
/*---------------------------------------------------------------------------*/
/* make correlation matrix by transforming symmetric positive definit matrix */
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

static void _unur_norta_debug_sigma_y( const struct unur_gen *gen, 
				       const double *sigma_y, 
				       const char *comment );
/*---------------------------------------------------------------------------*/
/* print sigma_y of corresponding normal distribution.                       */
/*---------------------------------------------------------------------------*/

static void _unur_norta_debug_eigensystem( const struct unur_gen *gen,
					   const double *eigenvalues,
					   const double *eigenvectors );
/*---------------------------------------------------------------------------*/
/* print eigensystem of sigma_y.                                             */
/*---------------------------------------------------------------------------*/

static void _unur_norta_debug_nmgenerator( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print genid of multinormal generator.                                     */
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

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_norta_debug_init(gen);
#endif

  /* compute parameters for NORTU (normal to uniform) */
  if (_unur_norta_nortu_setup(gen) != UNUR_SUCCESS) {
    free(par); _unur_norta_free(gen); return NULL;
  }

  /* distribution object for standard normal distribution */
  GEN.normaldistr = unur_distr_normal(NULL,0);




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

  /* allocate array for auxiliary copula */
  GEN.copula = _unur_xmalloc(sizeof(double)*GEN.dim);

  /* initialize pointer */
  NORMAL = NULL;
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
#define idx(a,b) ((a)*dim+(b))

  int dim = GEN.dim;    /* dimension of distribution */
  double *sigma_y;      /* correlation matrix for corresponding normal distr. */
  double *eigenvalues;  /* eigenvalues of sigma_y */
  double *eigenvectors; /* eigenvectors of sigma_y */
  double eigenvalues_positive; /* boolean indicating whether all eigenvalues are 
				  strictly positive */
  struct unur_distr *mn_distr; /* multinormal distribution */ 
  struct unur_gen   *mn_gen;   /* generator for multinormal distribution */ 
  int i,j;

  /* setup correlation matrix for corresponding normal distribution */
  sigma_y = _unur_xmalloc(dim * dim * sizeof(double));
  for(i=0; i<dim; i++) {
    /* left lower part: make matrix symmetric */
    for(j=0; j<i; j++)
      sigma_y[idx(i,j)] = sigma_y[idx(j,i)];
    /*   diagonal */
    sigma_y[idx(i,i)] = 1.;
    /* right upper part */
    for(j=i+1; j<dim; j++)
      sigma_y[idx(i,j)] = 2.*sin(DISTR.rankcorr[idx(i,j)]*(M_PI/6.));  
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
    _unur_norta_debug_sigma_y( gen, sigma_y, "NORTU setup:" );
#endif

  /* compute eigenvectors and eigenvalues of sigma_y */
  eigenvalues = _unur_xmalloc(dim * sizeof(double));
  eigenvectors = _unur_xmalloc(dim * dim * sizeof(double));
  if (_unur_matrix_eigensystem(dim, sigma_y, eigenvalues, eigenvectors) != UNUR_SUCCESS) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"cannot compute eigenvalues for given sigma_y");
    free(sigma_y); free(eigenvalues); free(eigenvectors);
    return UNUR_ERR_GEN_DATA;
  }

#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
    _unur_norta_debug_eigensystem( gen, eigenvalues, eigenvectors );
#endif

  /* check if all eigenvalues are positive */
  /* otherwise set to small values close to 0 */
  eigenvalues_positive = TRUE;
  for(i=0; i<dim; i++)
    if(eigenvalues[i] < UNUR_NORTA_MIN_EIGENVALUE) {
      eigenvalues[i] = UNUR_NORTA_MIN_EIGENVALUE;
      eigenvalues_positive = FALSE;
    }

  /* make corrected correlation matrix */
  if (!eigenvalues_positive) {
    _unur_matrix_transform_diagonal(dim,eigenvectors,eigenvalues,sigma_y);
    _unur_norta_make_correlationmatrix(dim,sigma_y);
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,
		  "sigma_y not positive definite -> corrected matrix");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
      _unur_norta_debug_sigma_y( gen, sigma_y, "\tEigenvalues < 0 --> correction required" );
#endif
  }

  /* clear working arrays */
  free(eigenvalues);
  free(eigenvectors);

  /* make generator for multinormal distribution */
  mn_distr = unur_distr_multinormal(dim, NULL, sigma_y);
  mn_gen = NULL;
  if (mn_distr) {
    mn_gen = unur_init(unur_vmt_new(mn_distr));
    _unur_distr_free(mn_distr);
  }
  if (mn_gen == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"(corrected) sigma_y not positive definit");
    free(sigma_y);
    return UNUR_ERR_GEN_DATA;
  }
  NORMAL = mn_gen;
  /* need same uniform random number generator as NORTA generator */
  NORMAL->urng = gen->urng;
  /* copy debugging flags */
  NORMAL->debug = gen->debug;

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
      _unur_norta_debug_nmgenerator( gen );
#endif

  /* clear working arrays */
  free(sigma_y);

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_norta_nortu_setup() */

/*---------------------------------------------------------------------------*/

int
_unur_norta_make_correlationmatrix( int dim, double *M)
     /*----------------------------------------------------------------------*/
     /* make correlation matrix by transforming symmetric positive definit   */
     /* matrix.                                                              */
     /*                                                                      */
     /* There is no checking whether M fulfills the conditions!              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim ... dimension of matrix M                                      */
     /*   M   ... symmetric square matrix                                    */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  int i,j;

  /* diagonal is used to store the roots of the diagonal elements */
  for (i=0; i<dim; i++)
    M[idx(i,i)] = sqrt(M[idx(i,i)]);

  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++)
      if(i!=j) M[idx(i,j)] /= M[idx(i,i)] * M[idx(j,j)];

  /* the diagonal elements are set to 1. */
  for (i=0; i<dim; i++) 
    M[idx(i,i)] = 1.;

  return UNUR_SUCCESS;
#undef idx
} /* end of _unur_norta_make_correlationmatrix() */

/*---------------------------------------------------------------------------*/

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
#define idx(a,b) ((a)*GEN.dim+(b))
  int j;
  double *u;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  /* pointer to auxiliary array of uniforms */
  u = GEN.copula;

  /* sample from multinormal distribution */
  _unur_sample_vec(NORMAL,u);

  /* make copula */
  for (j=0; j<GEN.dim; j++)
    u[j] = unur_distr_cont_eval_cdf( u[j], GEN.normaldistr );

  if (gen->distr->id == UNUR_DISTR_COPULA) {
    /* we want to have a normal copula --> just copy data */
    for (j=0; j<GEN.dim; j++) vec[j] = u[j];
    return;
  }



  fprintf(stderr,"junk\n");
  
  
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

/*   fprintf(log,"%s: data for corresponding normal distribution\n",gen->genid); */
  

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

void
_unur_norta_debug_sigma_y( const struct unur_gen *gen, 
			   const double *sigma_y, 
			   const char *comment )
     /*----------------------------------------------------------------------*/
     /* print sigma_y of corresponding normal distribution.                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   sigma_y ... pointer to correlation matrix                          */
     /*   comment ... additional string printed                              */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: %s\n",gen->genid,comment);
  fprintf(log,"%s:\n",gen->genid);
  _unur_matrix_print_matrix( GEN.dim, sigma_y, "\tsigma_y =", 
			     log, gen->genid, "\t   ");

} /* end of _unur_norta_debug_sigma_y() */

/*---------------------------------------------------------------------------*/

void
_unur_norta_debug_eigensystem( const struct unur_gen *gen,
			       const double *eigenvalues,
			       const double *eigenvectors )
     /*----------------------------------------------------------------------*/
     /* print eigensystem of sigma_y.                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*   eigenvalues  ... eigenvalues                                       */
     /*   eigenvectors ... eigenvalues                                       */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  log = unur_get_stream();

  _unur_matrix_print_vector( GEN.dim, eigenvalues, 
			     "\teigenvalues of sigma_y =", 
			     log, gen->genid, "\t   ");
  _unur_matrix_print_matrix( GEN.dim, eigenvectors, 
			     "\teigenvectors of sigma_y [rows] =", 
			     log, gen->genid, "\t   ");

} /* end of _unur_norta_debug_eigensystem() */

/*---------------------------------------------------------------------------*/

void
_unur_norta_debug_nmgenerator( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print genid of multinormal generator.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: generator for multinormal auxiliary distribution = %s\n", gen->genid,
	  NORMAL->genid );
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_norta_debug_nmgenerator() */

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
/* #define idx(a,b) ((a)*dim+(b)) */

/*   double *sigma_y;     /\* covariance matrix for ... *\/ */

  int i,j,notposdef,notpos;
  double *helpm,*eval,*evec;
  struct nortu_gen *gen;
  
/*   int dim = GEN.dim; */

  /* allocate working arrays */
/*   sigma_y = _unur_xmalloc(dim*dim*sizeof(double)); */


/*   gen = malloc(sizeof(struct nortu_gen)); */
/*   gen->dim = dim; */
  gen->cholesky = malloc(sizeof(double)*dim*dim);
/*   gen->distrN01 = unur_distr_normal(NULL,0);/\* to use the CDF in nortu_sample*\/ */
/*   gen->genN01 = unur_str2gen("normal()");/\*to generate N(0,1) variates in nortu_sample*\/ */
  helpm = malloc(sizeof(double)*dim*dim);

  //  printf("begin of nortusetup correlation matrix:\n");
  //  printsqmatrix(cor,dim);

  /* setup correlation matrix for ... */
/*   for(i=0; i<dim; i++) { */
/*     for(j=0; j<i; j++) */
/*       /\* left lower part: make matrix symmetric *\/ */
/* 	 sigma_y[idx(i,j)] = sigma_y[idx(j,i)]; */
/*     /\*   diagonal *\/ */
/*     sigma_y[idx(i,i)] = 1.; */
/*     for(j=i+1;j<dim;j++) */
/*       /\* right upper part *\/ */
/*       sigma_y[idx(i,j)] = 2.*sin(GEN.corr[idx(i,j)]*(M_PI/6.));   */
/*   } */


  //  printf("nortusetup sigma_y:\n");
  //  printsqmatrix(sigma_y,dim);

/*   gsl_matrix_view m */
/*     = gsl_matrix_view_array (sigma_y, dim, dim); */
/*   gsl_vector *evalgsl = gsl_vector_alloc (dim); */
/*   gsl_matrix *evecgsl = gsl_matrix_alloc (dim, dim); */
/*   gsl_eigen_symmv_workspace * w = */
/*     gsl_eigen_symmv_alloc (dim); */
/*   gsl_eigen_symmv (&m.matrix, evalgsl, evecgsl, w); */
/*   gsl_eigen_symmv_free (w); */
/*   eval = evalgsl->data; */
/*   evec = evecgsl->data; */

  //  printf("the eigenvalues\n");
  //  printvector(eval,dim);

  /*  printf("the eigenvectors\n");
  printsqmatrix(evec,dim);
  */
  /* the vector eval is now containing the eigenvalues
     and evec the corresponding eigenvectors */

  /*as gsl_eigen_symmv is destroying the diagonal and the lower triangular
    we have to restore the full matrix sigma_y  */
/*   for(i=0;i<dim;i++){ */
/*     for(j=0;j<i;j++) */
/*       sigma_y[idx(i,j)] = sigma_y[idx(j,i)]; */
/*     sigma_y[idx(i,i)] = 1.; */
/*   } */


  /* check if all eval are positive */
/*   notpos = 0; */
/*   for(i=0; i<dim; i++) */
/*     if(eval[i]<=0.){ /\*oder <= eps ?????*\/ */
/*       eval[i]=eps; */
/*       notpos++; */
/*     } */

/*   if(notpos>0){ */
/*     mmult_matr_diagonal(evec,eval,dim,helpm); */
/*     mmult_matr_matrt(helpm,evec,dim,sigma_y); */
/*     makecorrelationmatrix(sigma_y,dim); */
/*     //    printf("as sigma_y not positive definite EC corrected correlation matrix\n"); */
/*     //printsqmatrix(sigma_y,dim); */
/*   } */


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
/* #undef idx */
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
