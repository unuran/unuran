/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mcorr.c                                                      *
 *                                                                           *
 *   TYPE:      random matrix                                                *
 *   METHOD:    Matrix -- COORelation matrix                                 *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      random correlation matrix.                                           *
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
 *   [1] Devroye, L. (1986): Non-Uniform Random Variate Generation,          *
 *       New-York, Sect.6.1, p.605.                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Generate matrix H where all rows are independent and uniformly          *
 *   distributed on the sphere and return HH'.                               *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/matr.h>
#include <distributions/unur_distributions.h>
#include <uniform/urng.h>
#include <utils/unur_fp_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "arou.h"
#include "mcorr.h"

#include <utils/matrix_source.h>


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

#define MCORR_SET_EIGENVALUES  0x001u   /* set eigenvalues of corr-matrix    */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MCORR"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcorr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcorr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static void _unur_mcorr_sample_matr( struct unur_gen *gen, double *mat );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_mcorr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_mcorr_eigen (int dim, double *values, double *M );
/*---------------------------------------------------------------------------*/
/* Calculates a random correlation matrix M with given eigenvalues           */
/* using the Marsaglia-Olkin method                                          */
/*---------------------------------------------------------------------------*/


#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mcorr_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.matr      /* data for distribution object      */

#define PAR       par->data.mcorr       /* data for parameter object         */
#define GEN       gen->data.mcorr       /* data for generator object         */
#define DISTR     gen->distr->data.matr /* data for distribution in generator object */

#define SAMPLE    gen->sample.matr      /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/
#define NORMAL  gen->gen_aux        /* pointer to normal variate generator   */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_mcorr_new( const struct unur_distr *distr )
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
  if ( !(distr->type == UNUR_DISTR_MATR &&
	 distr->id == UNUR_DISTR_MCORRELATION) ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_MATR,NULL);

  /* allocate structure */
  par = _unur_xmalloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_MCORR_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  par->method   = UNUR_METH_MCORR;    /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* number of rows and columns (dimension of distribution). */
  /* do not confuse with distr->dim which is the size of     */
  /* the array that stores the matrix.                       */
  PAR.dim = distr->data.matr.n_rows;

  PAR.eigenvalues = NULL; /* (optional) eigenvalues of correlation matrix */

  /* routine for starting generator */
  par->init = _unur_mcorr_init;

  return par;

} /* end of unur_mcorr_new() */

/*****************************************************************************/

struct unur_gen *
_unur_mcorr_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_MCORR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MCORR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mcorr_create(par);
  if (!gen) { free(par); return NULL; }

  /* we need a generator for standard normal distributons */
  if (NORMAL==NULL) {
    struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
    struct unur_par   *normalpar = unur_arou_new( normaldistr );
    unur_arou_set_usedars( normalpar, TRUE );
    NORMAL = unur_init( normalpar );
    _unur_distr_free( normaldistr );
    if (NORMAL == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"Cannot create aux Gaussian generator");
      _unur_free(gen); free (par);
      return NULL;
    }
    /* need same uniform random number generator and debugging flags */
    NORMAL->urng = gen->urng;
    NORMAL->debug = gen->debug;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_mcorr_debug_init(gen);
#endif

  /* free parameters */
  if (PAR.eigenvalues) free(PAR.eigenvalues);
  free(par);

  /* o.k. */
  return gen;

} /* end of _unur_mcorr_init() */

/*---------------------------------------------------------------------------*/

int
unur_mcorr_set_eigenvalues( UNUR_PAR *par, double *eigenvalues  )
     /*----------------------------------------------------------------------*/
     /* sets the (optional) eigenvalues of the correlation matrix            */
     /*----------------------------------------------------------------------*/
{
  int i;
  double sum_eigenvalues = 0;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MCORR );
  CHECK_NULL( eigenvalues, UNUR_ERR_NULL );

  /* allocate memory if needed */
  if (PAR.eigenvalues==NULL)
    PAR.eigenvalues=_unur_xmalloc(PAR.dim * sizeof(double));

  /* set given eigenvalues */
  for (i=0; i<PAR.dim; i++) {
    PAR.eigenvalues[i] = eigenvalues[i];
    
    /* check for negative eigenvalues */
    if (PAR.eigenvalues[i]<0) {
      PAR.eigenvalues[i] = -PAR.eigenvalues[i];
      _unur_warning("MCORR", UNUR_ERR_GENERIC,"negative eigenvalue -> positive");
    }
    sum_eigenvalues += PAR.eigenvalues[i];
  }
  
  /* check if all eigenvalues = 0 ? */
  if (sum_eigenvalues==0) {
    _unur_error("MCORR", UNUR_ERR_GENERIC,"sum(eigenvalues)=0");
    return UNUR_ERR_GENERIC;
  }

  /* scaling values */
  if (!_unur_FP_same(sum_eigenvalues, (double) PAR.dim)) {
    _unur_warning("MCORR", UNUR_ERR_GENERIC,"scaling sum(eigenvalues) -> dim");
    for (i=0; i<PAR.dim; i++) {
      PAR.eigenvalues[i] = PAR.dim * PAR.eigenvalues[i]/sum_eigenvalues;
    }
  }

  /* changelog */
  par->set |= MCORR_SET_EIGENVALUES;

  return UNUR_SUCCESS;
}


/*---------------------------------------------------------------------------*/


static struct unur_gen *
_unur_mcorr_create( struct unur_par *par )
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
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MCORR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_MCORR_GEN);

  /* number of rows and columns (dimension of distribution). */
  /* do not confuse with distr->dim which is the size of     */
  /* the array that stores the matrix.                       */
  GEN.dim = DISTR.n_rows;

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_mcorr_sample_matr;
  gen->destroy = _unur_mcorr_free;
  gen->clone = _unur_mcorr_clone;

  /* allocate working array */
  GEN.H = _unur_xmalloc(GEN.dim * GEN.dim * sizeof(double));

  /* copy optional eigenvalues of the correlation matrix */
  GEN.eigenvalues = NULL;
  if (gen->set && MCORR_SET_EIGENVALUES) {
    GEN.eigenvalues = _unur_xmalloc(GEN.dim * sizeof(double));
    for (i=0; i<GEN.dim; i++) GEN.eigenvalues[i]=PAR.eigenvalues[i];
  }

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_mcorr_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mcorr_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.mcorr

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MCORR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate new working array */
  CLONE.H = _unur_xmalloc(GEN.dim * GEN.dim * sizeof(double));

  return clone;

#undef CLONE
} /* end of _unur_mcorr_clone() */

/*****************************************************************************/

void
_unur_mcorr_sample_matr( struct unur_gen *gen, double *mat )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   mat ... random matrix (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*(GEN.dim)+(b))
  int i,j,k;
  double sum, norm, x;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_MCORR_GEN,RETURN_VOID);
  CHECK_NULL(mat,RETURN_VOID);

  /* Marsaglia Olkin method with given eigenvalues */
  if (gen->set && MCORR_SET_EIGENVALUES) {
    _unur_mcorr_eigen(GEN.dim, GEN.eigenvalues, mat);
    return; 
  }

  /* eigenvalues have not been provided */
  /* generate rows vectors of matrix H uniformly distributed in the unit sphere */
  for (i=0; i<GEN.dim; i++) {
    sum=0.;
    for (j=0; j<GEN.dim; j++) {
      x = _unur_sample_cont(NORMAL);
      GEN.H[idx(i,j)] = x;
      sum += x * x;
    }
    norm = sqrt(sum);
    for (j=0; j<GEN.dim; j++) GEN.H[idx(i,j)] /= norm;
  }

  /* Compute HH' */
  for (i=0; i<GEN.dim; i++)
    for (j=0; j<GEN.dim; j++) {
      if (j<i)
	mat[idx(i,j)] = mat[idx(j,i)];
      else if(j==i)
	mat[idx(i,j)] = 1.;
      else {
	sum=0.;
	for (k=0; k<GEN.dim; k++)
	  sum += GEN.H[idx(i,k)]*GEN.H[idx(j,k)];
	mat[idx(i,j)] = sum;
      }
    }

#undef idx
} /* end of _unur_mcorr_sample_matr() */

/*****************************************************************************/

void
_unur_mcorr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_MCORR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MCORR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN.eigenvalues) free(GEN.eigenvalues);
  if (GEN.H)           free(GEN.H);

  _unur_generic_free(gen);

} /* end of _unur_mcorr_free() */

/*--------------------------------------------------------------------------*/

int _unur_mcorr_eigen (int dim, double *values, double *M )
     /* Calculates a random correlation matrix M with given eigenvalues     */
     /* using the Marsaglia-Olkin method                                    */
{
#define idx(a,b) ((a)*dim+(b))
  int i,j,k;
  double *E, *P;
  double *x, *y, *z, *w, *r; /* misc vectors used in the marsaglia-olkin method */
  double a, b, c, e, e2;
  int s; /* random sign +-1 */

  /* check parameters */
  CHECK_NULL(values, UNUR_ERR_NULL);
  CHECK_NULL(M, UNUR_ERR_NULL);
  if (dim<1) {
    _unur_error("MCORR",UNUR_ERR_GENERIC,"dimension < 1");
    return UNUR_ERR_GENERIC;
  }

  /* initialization steps */
  x=_unur_xmalloc(dim*sizeof(double));
  y=_unur_xmalloc(dim*sizeof(double));
  z=_unur_xmalloc(dim*sizeof(double));
  w=_unur_xmalloc(dim*sizeof(double));
  r=_unur_xmalloc(dim*sizeof(double));
  E=_unur_xmalloc(dim*dim*sizeof(double));
  P=_unur_xmalloc(dim*dim*sizeof(double));

  /* initially E is an identity matrix */
  for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    E[idx(i,j)] = (i==j) ? 1: 0;
  }}


  for (k=0; k<dim-1; k++) {
    /* w is a random vector */
    for (i=0; i<dim; i++) w[i] = unur_urng_sample(NULL);

    /* x = E*w */
    for (i=0; i<dim; i++) {
      x[i]=0;
      for (j=0; j<dim; j++) {
        x[i] += E[idx(i,j)]*w[j];
      }
    }

    /* a=sum{(1-lambda_i)*x_i*x_i} */
    a=0;
    for (i=0; i<dim; i++)
      a += (1-values[i])*x[i]*x[i];

    /* check if all eigenvalues are ~1 */
    if (fabs(a)<DBL_EPSILON) {
      /* return identity matrix */
      for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) {
        M[idx(i,j)] = (i==j) ? 1: 0;
      }}
      _unur_warning("MCORR", UNUR_ERR_GENERIC,"all eigenvalues are ~1 -> identity matrix");
      
      return UNUR_SUCCESS;
    }

    do {

      /* z is a random vector */
      for (i=0; i<dim; i++) z[i] = unur_urng_sample(NULL);

      /* y = E*z */
      for (i=0; i<dim; i++) {
        y[i]=0;
        for (j=0; j<dim; j++) {
          y[i] += E[idx(i,j)]*z[j];
        }
      }

      /* b=sum{(1-lambda_i)*x_i*y_i} */
      /* c=sum{(1-lambda_i)*y_i*y_i} */
      b=0; c=0;
      for (i=0; i<dim; i++) {
        b += (1-values[i])*x[i]*y[i];
        c += (1-values[i])*y[i]*y[i];
      }

      /* e^2 = b^2 - a*c */
      e2 = b*b - a*c;

    } while (e2<0);

    e=sqrt(e2);


    /* random sign */
    s = (unur_urng_sample(NULL)>.5) ? 1: -1 ;

    /* r=x*(b+s*e)/a - y */
    for (i=0; i<dim; i++) r[i] = x[i]*(b+s*e)/a - y[i];

    /* another random sign */
    s = (unur_urng_sample(NULL)>.5) ? 1: -1 ;

    /* pk=s*r/norm(r) */
    _unur_vector_normalize(dim, r);
    for (i=0; i<dim; i++) P[idx(k,i)] = s * r[i];

    /* E = E - r r^T */
    for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) {
        E[idx(i,j)] -= r[i]*r[j];
      }
    }

  } /* next k */

  /* w is a random vector */
  for (i=0; i<dim; i++) w[i] = unur_urng_sample(NULL);

  /* x = E*w */
  for (i=0; i<dim; i++) {
    x[i]=0;
    for (j=0; j<dim; j++) {
      x[i] += E[idx(i,j)]*w[j];
    }
  }

  _unur_vector_normalize(dim, x);

  /* last row of the orthogonal matrix P */
  for (i=0; i<dim; i++) {
    P[idx(dim-1,i)] = x[i];
  }

  /* M = P L P^T, where L diagonal containing the eigenvalues */
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      M[idx(i,j)] = 0;
      for (k=0; k<dim; k++) {
        M[idx(i,j)] += P[idx(i,k)] * values[k] * P[idx(j,k)];
      }
    }
  }

  free(E); free(P);
  free(x); free(y); free(z); free(w); free(r);

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_mcorr_eigen() */

/*--------------------------------------------------------------------------*/





/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_mcorr_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MCORR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = random matrix\n",gen->genid);
  fprintf(log,"%s: method  = MCORR (Matrix - CORRELATION matrix)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_matr_debug( gen->distr, gen->genid );

  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = _unur_mcorr_sample_matr()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_mcorr_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
