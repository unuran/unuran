/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      varou.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    adaptive ratio-of-uniforms method                            *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      ...                                                                  *
 *      Produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      PDF, dPDF                                                            *
 *   OPTIONAL:                                                               *
 *      ...                                                                  *
 *                                                                           *
 *****************************************************************************
 *   (c) 2000 Wolfgang Hoermann and Josef Leydold                            *
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
 *   [1] Wakefield J.C., Gelfand A.E., Smith A.F.M.                          *
 *       Efficient generation of random variates via the ratio-of-uniforms   *
 *       method.                                                             *
 *       Statistics and Computing (1991) 1, pp (129-133)                     *
 *                                                                           *
 *   [2] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <utils/fmax_source.h>
#include <utils/hooke_source.h> 
#include <utils/matrix_source.h>
#include <utils/vector_source.h>
#include <utils/unur_fp_source.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "varou.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define VAROU_VARFLAG_VERIFY   0x002u   /* run verify mode                   */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define VAROU_SET_SOMETHING    0x001u     /* set something                   */

/*---------------------------------------------------------------------------*/

#define GENTYPE "VAROU"         /* type of generator                         */

/*---------------------------------------------------------------------------*/

/* convergence parameters for the hooke optimization algorithm */

#define VAROU_HOOKE_RHO 0.5 
#define VAROU_HOOKE_EPSILON 1e-7 
#define VAROU_HOOKE_MAXITER 10000

/* scaling factor of computed minimum boundary rectangle */
/* after we have computed the boundary rectangle (0, vmax)x(umin[d], umax[d])*/
/* we scale the obtained boundaries with this factor, i.e. :                 */
/* vmax = vmax * ( 1+ VAROU_RECT_SCALING)                                    */
/* umin[d] = umin[d] - (umax[d]-umin[d])*VAROU_RECT_SCALING/2.               */
/* umax[d] = umax[d] + (umax[d]-umin[d])*VAROU_RECT_SCALING/2.               */
#define VAROU_RECT_SCALING 1e-4

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_varou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_varou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static void  _unur_varou_sample_cvec( struct unur_gen *gen, double *vec );
static void  _unur_varou_sample_check( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static int _unur_varou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static void _unur_varou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/


#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_varou_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif


/*---------------------------------------------------------------------------*/
/*   Auxilliary Routines                                                     */
/*---------------------------------------------------------------------------*/

double *_unur_varou_sample_simplex(struct unur_gen *gen, double *vertices);
/*---------------------------------------------------------------------------*/
/* Generate random vector uniformly distributed in simplex                   */
/*---------------------------------------------------------------------------*/


double _unur_varou_aux_vmax(double *x, void *p );
double _unur_varou_aux_umin(double *x, void *p );
double _unur_varou_aux_umax(double *x, void *p );
/*---------------------------------------------------------------------------*/
/* Auxiliary functions used in the computation of the bounding rectangle     */
/*---------------------------------------------------------------------------*/

static double _unur_varou_f( struct unur_gen *gen, double *u, double v);
/*---------------------------------------------------------------------------*/
/* return calculated value of PDF(u_1/v, u_2/v, ..., u_dim/v)                */
/*---------------------------------------------------------------------------*/

static double _unur_varou_F( struct unur_gen *gen, double *u, double v);
/*---------------------------------------------------------------------------*/
/* return calculated value of v^(1+dim) - PDF(u_1/v, u_2/v, ..., u_dim/v)    */
/*---------------------------------------------------------------------------*/

static double *_unur_varou_dF( struct unur_gen *gen, double *u, double v);
/*---------------------------------------------------------------------------*/
/* return calculated value the gradient of F() at the point (u,v)            */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/

/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       par->data.varou        /* data for parameter object        */
#define GEN       gen->data.varou        /* data for generator object        */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */     
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_varou_new( const struct unur_distr *distr )
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

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }

  /* allocate structure */
  par = _unur_xmalloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_VAROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* copy number of dimensions from the distribution object */
  PAR.dim = distr->dim;

  /* set default values */
  par->method   = UNUR_METH_VAROU;    /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_varou_init;

  return par;

} /* end of unur_varou_new() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_varou_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VAROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | VAROU_VARFLAG_VERIFY) : 
                            (par->variant & (~VAROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_varou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_varou_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VAROU, UNUR_ERR_GEN_INVALID );

  if (verify) {
    /* turn verify bounding rectangle on */
    gen->variant |= VAROU_VARFLAG_VERIFY;
    SAMPLE = _unur_varou_sample_check;
  }
  else {
    /* turn verify bounding rectangle off */
    gen->variant &= ~VAROU_VARFLAG_VERIFY;
    SAMPLE = _unur_varou_sample_cvec;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_varou_chg_verify() */

/*****************************************************************************/

struct unur_gen *
_unur_varou_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params ... pointer to paramters for building generator object      */
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
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_VAROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VAROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_varou_create(par);
  if (!gen) { 
    free(par); 
    return NULL; 
  }

  /* compute bounding rectangle */
  if (_unur_varou_rectangle(gen)!=UNUR_SUCCESS) {
    free(par); _unur_varou_free(gen);
    return NULL;
  }


  /* compute initial triangulation ... etc  */
  /* ...  */


#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_varou_debug_init(gen);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of _unur_varou_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_varou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VAROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_VAROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & VAROU_VARFLAG_VERIFY) ? _unur_varou_sample_check : 
                                                   _unur_varou_sample_cvec;

  gen->destroy = _unur_varou_free;
  gen->clone = _unur_varou_clone;

  /* allocate memory for u-boundary arrays */
  GEN.umin = _unur_xmalloc( PAR.dim * sizeof(double)); /* bounding rectangle */
  GEN.umax = _unur_xmalloc( PAR.dim * sizeof(double)); /* bounding rectangle */

  /* copy parameters into generator object */
  GEN.dim   = PAR.dim;              /* dimension */
  
  /* get center of the distribution */
  GEN.center = unur_distr_cvec_get_center(gen->distr);

  /* initialize parameters */

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_varou_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_varou_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.vnrou

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VAROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate memory for u-arrays and center */
  CLONE.umin = _unur_xmalloc( GEN.dim * sizeof(double));
  CLONE.umax = _unur_xmalloc( GEN.dim * sizeof(double));
  
  /* copy parameters into clone object */
  memcpy(CLONE.umin, GEN.umin, GEN.dim * sizeof(double));
  memcpy(CLONE.umax, GEN.umax, GEN.dim * sizeof(double));

  /* copy data */
  CLONE.center = unur_distr_cvec_get_center(clone->distr);
  
  return clone;

#undef CLONE
} /* end of _unur_varou_clone() */

/*****************************************************************************/

void
_unur_varou_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{ 
  int d, dim; /* index used in dimension loops (0 <= d < dim) */

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  
  COOKIE_CHECK(gen,CK_VAROU_GEN,RETURN_VOID); 

  dim = GEN.dim;
 
  /* TODO : CHANGE THIS !!!!!!! These are only dummy random values */  
    for (d=0; d<dim; d++) {
      vec[d] = _unur_call_urng(gen->urng) ;
    }
    
  return; 
   
} /* end of _unur_varou_sample() */

/*---------------------------------------------------------------------------*/

void
_unur_varou_sample_check( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random sample vector (return)                              */
     /*----------------------------------------------------------------------*/
{ 
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  
  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  
  COOKIE_CHECK(gen,CK_VAROU_GEN,RETURN_VOID); 

  dim = GEN.dim;
  
  /* TODO : CHANGE THIS !!!!!!! These are only dummy random values */  
    for (d=0; d<dim; d++) {
      vec[d] = _unur_call_urng(gen->urng) ;
    }
 
  return;

} /* end of _unur_varou_sample_check() */

/*****************************************************************************/

void
_unur_varou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_VAROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VAROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN.umin) free(GEN.umin); 
  if (GEN.umax) free(GEN.umax);

  _unur_generic_free(gen);

} /* end of _unur_varou_free() */

/*****************************************************************************/

double *
_unur_varou_sample_simplex( struct unur_gen *gen, double *vertices)
    /*----------------------------------------------------------------------*/
    /* Generate random vector uniformply distributed in simplex             */
    /* having (dim+1) vertices (v_0, ... , v_dim)                           */
    /* The coordinates of the (dim+1) vertices are stored sequentially in   */
    /* the vertices[] array of size (dim+1)*dim                             */
    /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  double *X; /* random sample vector */
  double  U; /* uniform variate */
  double *E; /* exponential variates */
  double *S; /* uniform spacings */
  double  E_sum; /* sum of all e[] */
  long i,j; /* index variables used in for loops */
  long dim;
  
  dim = GEN.dim; 
  
  X=_unur_xmalloc(dim*sizeof(double));
  S=_unur_xmalloc((dim+1)*sizeof(double));
  E=_unur_xmalloc((dim+1)*sizeof(double));
  
  /* calculating exponential variates */
  E_sum = 0; 
  for (i=0; i<=dim; i++) {
     /* sample from U(0,1) */
     while ( (U = _unur_call_urng(gen->urng)) == 0.);
     /* sample from exponential distribution */
     E[i] = -log(1.-U);
     E_sum += E[i];	 
  }
  
  /* calculating uniform spacing */
  for (i=0; i<=dim; i++) {
     S[i]=E[i]/E_sum;
  }

  /* calculating uniform vector in simplex */
  for (j=0; j<dim; j++) {
    X[j] = 0;
    for (i=0; i<=dim; i++) {
      X[j] += S[i] * vertices[idx(i,j)];
    }
  }

  free(S); free(E);

  return X;

#undef idx
} /* end of _unur_varou_sample_simplex() */

/*****************************************************************************/



/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/
 
double
_unur_varou_aux_vmax(double *x, void *p )
    /*----------------------------------------------------------------------*/
    /* Auxiliary function used in the computation of the bounding rectangle */
    /*----------------------------------------------------------------------*/
{

    struct unur_gen *gen;
    gen = p; /* typecast from void* to unur_gen* */

    return -pow( _unur_cvec_PDF((x),(gen->distr)) ,
                         1./(1.+ GEN.dim) );

} /* end of _unur_varou_aux_vmax() */
		   
/*---------------------------------------------------------------------------*/

double
_unur_varou_aux_umin(double *x, void *p)
    /*----------------------------------------------------------------------*/
    /* Auxiliary function used in the computation of the bounding rectangle */
    /*----------------------------------------------------------------------*/
{
    struct unur_gen *gen;

    gen = p; /* typecast from void* to unur_gen* */

    return (x[GEN.aux_dim] - GEN.center[GEN.aux_dim])
	              * pow( _unur_cvec_PDF((x),(gen->distr)),
		                1. / (1.+ GEN.dim) );

} /* end of _unur_varou_aux_umin() */
					      
/*---------------------------------------------------------------------------*/

double
_unur_varou_aux_umax(double *x, void *p)
    /*----------------------------------------------------------------------*/
    /* Auxiliary function used in the computation of the bounding rectangle */
    /*----------------------------------------------------------------------*/
{
    return (- _unur_varou_aux_umin(x,p)) ;

} /* end of _unur_varou_aux_umin() */

/*---------------------------------------------------------------------------*/


double _unur_varou_f( struct unur_gen *gen, double *u, double v) 
    /*                                                            */
    /* return calculated value of PDF(u_1/v, u_2/v, ..., u_dim/v) */
    /*                                                            */
    /* UNUR_INFINITY is returned when |v|<UNUR_EPSILON            */
    /*                                                            */
{
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double *uv; /* u/v vector     */
  double f;   /* function value */

  if (fabs(v) <= UNUR_EPSILON) return UNUR_INFINITY;
  
  dim = GEN.dim;

  uv = _unur_vector_new(dim);
  for (d=0; d<dim; d++) {
    uv[d] = u[d] / v;
  }

  f = PDF(uv);

  free(uv);
  
  return f;
}

/*---------------------------------------------------------------------------*/

double 
_unur_varou_F( struct unur_gen *gen, double *u, double v) 
    /*                                                                        */
    /* return calculated value of v^(1+dim) - PDF(u_1/v, u_2/v, ..., u_dim/v) */
    /*                                                                        */
    /* UNUR_INFINITY is returned when |v|<UNUR_EPSILON                        */
    /*                                                                        */
{
  double F;   /* function value */
 
  if (fabs(v) <= UNUR_EPSILON) return UNUR_INFINITY;
  
  F = pow(v, 1.+ GEN.dim ) - _unur_varou_f( gen, u, v);

  return F;
}

/*---------------------------------------------------------------------------*/

double * 
_unur_varou_dF( struct unur_gen *gen, double *u, double v )
    /*                                                                        */
    /* return the calculated value of the gradient of F() at the point (u,v)  */
    /*                                                                        */
    /* NULL is returned when |v|<UNUR_EPSILON                                 */
    /*                                                                        */
{
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double *uv; /* u/v vector     */
  double *dF;
  double *dPDF;

  if (fabs(v) <= UNUR_EPSILON) return NULL;

  dim = GEN.dim;
   
  /* reserving memory for arrays */ 
  uv = _unur_vector_new(dim);
  dPDF = _unur_vector_new(dim);
  dF = _unur_vector_new(dim+1);
  
  /* calculating the rations uv[d] = u[d]/v */
  for (d=0; d<dim; d++) {
    uv[d] = u[d] / v;
  }

  /* dPDF at the point uv[] */
  _unur_cvec_dPDF(dPDF, uv, gen->distr);

  /* gradient coordinates of F() */
  for (d=0; d<dim; d++) {
    dF[d] = - dPDF[d] / v;
  }
  /* and the last coordinate */
  dF[dim] = (1.+dim) * pow(v, dim);

  /* freeing allocated memory */
  _unur_vector_free(uv);
  _unur_vector_free(dPDF);
  
  return dF;
}

/*---------------------------------------------------------------------------*/

int
_unur_varou_rectangle( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 

  struct unur_funct_vgeneric faux; /* function to be minimized/maximized    */
  double *xstart, *xend, *xumin, *xumax; /* coordinate arrays used in maximum/minimum calculations */
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  int hooke_iters_vmax;  /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umin;  /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umax;  /* actual number of min/max iterations = return value of hooke()*/
  double scaled_epsilon; /* to be used in the hooke algorithm */

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_VAROU_GEN, UNUR_ERR_COOKIE );

  /* dimension of the distribution */
  dim = GEN.dim;

  /* allocate memory for the coordinate vectors */
  xstart = _unur_xmalloc(dim * sizeof(double));
  xend   = _unur_xmalloc(dim * sizeof(double));
  xumin  = _unur_xmalloc(dim * sizeof(double));
  xumax  = _unur_xmalloc(dim * sizeof(double));
  
  /* calculation of vmax */
    /* calculating vmax as maximum of sqrt(f(x)) in the domain */
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_varou_aux_vmax;
      faux.params = gen;

      if (gen->distr->set & UNUR_DISTR_SET_MODE)
          /* starting point at mode */
          memcpy(xstart, DISTR.mode, dim * sizeof(double)); 
      else
          /* starting point at center */
          memcpy(xstart, GEN.center, dim * sizeof(double)); 
     
      hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend, 
				      VAROU_HOOKE_RHO, VAROU_HOOKE_EPSILON, VAROU_HOOKE_MAXITER);

      GEN.vmax = -faux.f(xend, faux.params);
      
      if (hooke_iters_vmax >= VAROU_HOOKE_MAXITER) {
	 scaled_epsilon = VAROU_HOOKE_EPSILON * GEN.vmax;
	 if (scaled_epsilon>VAROU_HOOKE_EPSILON) scaled_epsilon=VAROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         memcpy(xstart, xend, dim * sizeof(double)); 
         hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend, 
                              VAROU_HOOKE_RHO, scaled_epsilon , VAROU_HOOKE_MAXITER);
         GEN.vmax = -faux.f(xend, faux.params);
         if (hooke_iters_vmax >= VAROU_HOOKE_MAXITER) {
           _unur_warning(gen->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (vmax)");  
         }
      }

      /* (u,v)-coordinates are (0,...,0,vmax) */


  /* calculation of umin and umax */
    
    for (d=0; d<dim; d++) {

      /* setting coordinate dimension to be used by the auxiliary functions */
      GEN.aux_dim  = d;
  
      if (gen->distr->set & UNUR_DISTR_SET_MODE)
          /* starting point at mode */
          memcpy(xstart, DISTR.mode, dim * sizeof(double)); 
      else
          /* starting point at center */
          memcpy(xstart, GEN.center, dim * sizeof(double)); 
      
      /*-----------------------------------------------------------------------------*/
      /* calculation for umin */
      
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_varou_aux_umin;
      faux.params = gen;

      hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend, 
                           VAROU_HOOKE_RHO, VAROU_HOOKE_EPSILON, VAROU_HOOKE_MAXITER);
      GEN.umin[d] = faux.f(xend, faux.params);
      
      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumin, xend, dim * sizeof(double)); 

      /*-----------------------------------------------------------------------------*/
      /* and now, an analogue calculation for umax */

      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_varou_aux_umax;
      faux.params = gen;

      hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend, 
                           VAROU_HOOKE_RHO, VAROU_HOOKE_EPSILON, VAROU_HOOKE_MAXITER);
      GEN.umax[d] = -faux.f(xend, faux.params);
      
      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumax, xend, dim * sizeof(double)); 

      /*-----------------------------------------------------------------------------*/
      /* checking if we need to recalculate umin */
      if (hooke_iters_umin >= VAROU_HOOKE_MAXITER) {
	 scaled_epsilon = VAROU_HOOKE_EPSILON * (GEN.umax[d]-GEN.umin[d]);
	 if (scaled_epsilon>VAROU_HOOKE_EPSILON) scaled_epsilon=VAROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         faux.f = (UNUR_FUNCT_VGENERIC*) _unur_varou_aux_umin;
         faux.params = gen;

         memcpy(xstart, xumin, dim * sizeof(double)); 
         hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend, 
                              VAROU_HOOKE_RHO, scaled_epsilon , VAROU_HOOKE_MAXITER);
         GEN.umin[d] = faux.f(xend, faux.params);
         if (hooke_iters_umin >= VAROU_HOOKE_MAXITER) {
           _unur_warning(gen->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umin)");  
         }

         /* storing actual endpoint */
         memcpy(xumin, xend, dim * sizeof(double));

      }

      /* checking if we need to recalculate umax */
      if (hooke_iters_umax >= VAROU_HOOKE_MAXITER) {
	 scaled_epsilon = VAROU_HOOKE_EPSILON * (GEN.umax[d]-GEN.umin[d]);
	 if (scaled_epsilon>VAROU_HOOKE_EPSILON) scaled_epsilon=VAROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         faux.f = (UNUR_FUNCT_VGENERIC*) _unur_varou_aux_umax;
         faux.params = gen;
        
	 memcpy(xstart, xumax, dim * sizeof(double)); 
         hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend, 
                              VAROU_HOOKE_RHO, scaled_epsilon , VAROU_HOOKE_MAXITER);
         GEN.umin[d] = faux.f(xend, faux.params);
         if (hooke_iters_umax >= VAROU_HOOKE_MAXITER) {
           _unur_warning(gen->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umax)");  
         }
         
	 /* storing actual endpoint */
         memcpy(xumax, xend, dim * sizeof(double));

         /* TODO: calculation of (u_i,v)-coordinates for u_i <> u_d */


      }
    
      /*-----------------------------------------------------------------------------*/
      /* additional scaling of boundary rectangle */   
      GEN.vmax = GEN.vmax * ( 1+ VAROU_RECT_SCALING);
      GEN.umin[d] = GEN.umin[d] - (GEN.umax[d]-GEN.umin[d])*VAROU_RECT_SCALING/2.;
      GEN.umax[d] = GEN.umax[d] + (GEN.umax[d]-GEN.umin[d])*VAROU_RECT_SCALING/2.;        
    
    }
 
  free(xstart); free(xend); free(xumin); free(xumax);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_varou_rectangle() */

/*---------------------------------------------------------------------------*/





/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_varou_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double vol;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VAROU_GEN,RETURN_VOID);

  log = unur_get_stream();
  dim = GEN.dim;

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = varou (adaptive ratio-of-uniforms)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  
  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_varou_sample",gen->genid);
  if (gen->variant & VAROU_VARFLAG_VERIFY) fprintf(log,"_check");
  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */


  /* print center */
  _unur_matrix_print_vector( GEN.dim, GEN.center, "center =", log, gen->genid, "\t   ");

  /* print bounding rectangle */
  fprintf(log,"%s: Rectangle:",gen->genid);
  fprintf(log,"\t[computed]\n");
 
  vol = GEN.vmax;
  fprintf(log,"%s:\tvmax = %g\n",gen->genid, GEN.vmax);
  for (d=0; d<dim; d++) {
    vol *= (GEN.umax[d]-GEN.umin[d]);
    fprintf(log,"%s:\tumin[%d],umax[%d] = (%g,%g)\n",gen->genid, 
	    d, d, GEN.umin[d], GEN.umax[d]);
  }
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s:\tvolume = %g\t(hat = %g)\n",gen->genid, vol, vol*(GEN.dim+1));
  fprintf(log,"%s:\n",gen->genid);


  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_varou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
