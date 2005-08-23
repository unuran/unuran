/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *                        -- EXPERIMENTAL CODE --                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl_fast_init.ch                                            *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a unimodal distribution                                 *
 *      produce random variate X with its density                            *
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
/* Constants                                                                 */

#define TABL_FAST_BLKSIZE  100
/* 
   Size block to allocated for the double array used by the fast table
   method.
*/

int unur_tabl_set_variant_fast( UNUR_PAR *parameters );
/* 
   Use a fast table version is used
   which is built only using the ``equal area rule'', without
   (derandomized) adaptive rejection sampling. (Other settings given
   by the user are ignored.) It uses a double array. Thus marginal
   generation times are very fast.  
   However, the areafraction is the only parameter when this variant
   is chosen which might result a suboptimal hat function. There is no
   adaptive algorithm to adjust the hat to the given PDF.
*/

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_tabl_gen { 
  double  Atotal;               /* total area below hat                      */
  double  Asqueeze;             /* area of squeeze polygon                   */

  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */

  double  *tabl;                /* pointer to array of data (for fast table method) */
  int     tabl_size;            /* size of data array                        */

  struct unur_tabl_interval **guide; /* pointer to guide table               */
  int     guide_size;           /* size of guide table                       */
  double  guide_factor;         /* relative size of guide table              */

  double  Umin, Umax;           /* bounds for iid random variable in respect to
				   the given (truncated) domain of the distr.*/

  struct unur_tabl_interval *iv;     /* pointer to linked list of intervals  */
  int     n_ivs;                /* number of intervals                       */
  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* limit for ratio r_n = A(squeeze) / A(hat) */

  double  darsfactor;           /* factor for (derandomized) ARS             */
};

/*---------------------------------------------------------------------------*/
/* The generator object for fast table method                                */

struct unur_tabl_fast_gen { 
  int     n_slopes;             /* number of slopes (given or computed)      */
  double  *slopes;              /* pointer to array of slopes                */

  double  *table;               /* pointer to array of data                  */
  int     n_table;              /* size of data array                        */

  double  bararea;              /* area of bars                              */
  double  n_ivs;                /* total number of intervals                 */  

  double  Atotal;               /* total area below hat                      */
  double  Asqueeze;             /* area of squeeze polygon                   */

  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */
};

/*---------------------------------------------------------------------------*/

static double _unur_tabl_fast_sample( struct unur_gen *gen );
static double _unur_tabl_fast_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_fast_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_fast_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_fast_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a single array for storing data and remove linked list.              */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_fast_debug_init_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_fast_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

#define fastGEN   ((struct unur_tabl_fast_gen*)gen->datap) /* data for fast table method */

/*****************************************************************************/
/**  Public                                                                 **/
/*****************************************************************************/

int
unur_tabl_set_variant_fast( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Switch between linked list and fast table                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* we use a bit in variant */
  par->variant = (par->variant & ~TABL_VARMASK_VARIANT) | TABL_VARIANT_FAST;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tabl_set_variant_fast() */

/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

int
_unur_tabl_fast_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a single array for storing data and remove linked list.         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define add_to_table(x) ( ((j<n_table)||(table=_unur_xrealloc(table,(n_table+=TABL_FAST_BLKSIZE)*sizeof(double)))) \
                            && (table[j++]=(x)))

  struct unur_tabl_interval *iv, *next;
  struct unur_tabl_gen *genslopes;
  double *table = NULL;
  int n_table = 0;
  double x, b, c;
  int i, j, k, n;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_GEN_INVALID );
  
  /* pointer to (old) data structure */
  /* has been used to store the slopes */
  genslopes = GEN;

  /* allocate memory for new data structure */
  gen->datap = _unur_xmalloc( sizeof(struct unur_tabl_fast_gen) );
  gen->s_datap = sizeof(struct unur_tabl_fast_gen);

  /* magic cookies */
  COOKIE_SET(gen,CK_TABL_fast_GEN);

  /* routines for cloning and destroying generator */
  gen->destroy = _unur_tabl_fast_free;
  gen->clone = _unur_tabl_fast_clone;

  /* set number of slopes and allocate memory */
  fastGEN->n_slopes = PAR->n_slopes;
  fastGEN->slopes = _unur_xmalloc( 2 * fastGEN->n_slopes * sizeof(double) );

  /* area for each bar */
  fastGEN->bararea = DISTR.area * PAR->area_fract;

  /* initialize double array */
  fastGEN->table = NULL;
  fastGEN->n_table = 0;

  /* counter for intervals */
  j = 0;
  k = 0;

  /* copy data about slopes into new data structure and  */
  /* remove entries from linked list.                    */
  for (iv = genslopes->iv, i=0; iv!=NULL; iv=next, i++) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_SHOULD_NOT_HAPPEN);
    if (i >= fastGEN->n_slopes) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_GENERIC;
    }

    /* copy data */
    (fastGEN->slopes)[2*i]   = iv->xmax;
    (fastGEN->slopes)[2*i+1] = iv->xmin;

    /* make intervals */
    n = j++;
    x = (fastGEN->slopes)[2*i];
    b = (fastGEN->slopes)[2*i+1];
    c = fastGEN->bararea * (x<b ? 1.:-1.); 
    if (x < b) { /* slope decreasing */
      while (_unur_FP_less(x,b)) { add_to_table(x); x += c/PDF(x); }
    }
    else { /* slope increasing */
      while (_unur_FP_greater(x,b)) { add_to_table(x); x += c/PDF(x); }
    }

    /* end of slope */
    if (_unur_isfinite(x)) {
      /* terminating interval */
      add_to_table(x);
      /* extra interval to store squeeze */
      x += c/iv->fmin;
      add_to_table(x);
    }
    else { 
      /* terminating interval has probability 0  */
      /* we have no squeeze                      */
      add_to_table(c*INFINITY);
    }

    /* store end of slope */
    add_to_table(b);

   /* store last regular interval of slope and update total number of intervals */
    table[n] = (double)(j-4);
    k += j-4 - n;

    /* free interval */
    next = iv->next;
    free(iv);
  }
  /* remove old data structure */
  free(genslopes);

  /* resize allocated memory block for table and store in data structure */
  fastGEN->n_table = j;
  fastGEN->table = _unur_xrealloc(table, j*sizeof(double));
  fastGEN->n_ivs = (double) k;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tabl_fast_debug_init_finished( gen );
#endif

  /* o.k. */
  return UNUR_SUCCESS;

#undef add_to_table
} /* end of _unur_tabl_fast_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tabl_fast_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_tabl_fast_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TABL_fast_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_tabl_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_tabl_fast_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_TABL ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TABL_fast_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tabl_fast_debug_free(gen);
#endif

  /* free table */
  if (fastGEN->slopes) free(fastGEN->slopes);
  if (fastGEN->table)  free(fastGEN->table);

  /* free other memory */
  _unur_generic_free(gen);

} /* end of _unur_tabl_fast_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double
_unur_tabl_fast_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use fast table)                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double *table;
  double r,s,w;
  double d,h,q;
  double fmin, fmax;
  double X,U,V;
  int j,m,n;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TABL_fast_GEN,INFINITY);

  table = fastGEN->table;
  r = fastGEN->bararea;
  s = fastGEN->n_ivs;
  
  do {
    /* w ~ U(0,s) */
    w = s * _unur_call_urng(gen->urng);  
    /* find slope (and interval) by sequential search */
    j = m = 1 + (int) w;
    n = 0;
    do {
      n += table[n];
      if (j <= n) break;
      n += 4, j+= 4;
    } while (1);
    /* compute width of interval and of "squeeze interval" */
    d = table[j+1]-table[j];     /* length of interval           = area / fmax */
    h = table[j+2]-table[j+1];   /* length of "squeeze interval" = area / fmin */
    U = w-m+1.;                  /* U ~ U(0,1) */
    q = h * U;                   /* q ~ U(0,h) */
    /* below or above squeeze */
    if (U <= d/h) {
      /* below squeeze --> immediate acceptance */
      X = table[j] + h * U;     /* h*U ~ U(0,d) --> X ~ U(a_j, a_{j+1}) */
      if (j<n) return X;    /* interval entirely inside slope --> accept */
      if ( (d>0 && X <= table[n+3]) || (d<0 && X >= table[n+3]) ) 
	/* accept only points in terminal interval which are inside of slope */
	return X;
    }
    else {
      /* reject from region below hat and squeeze */
      fmin = _unur_isfinite(h) ? r/fabs(h) : 0.; 
      fmax = r/fabs(d);
      q = fmin/fmax;

      /* compute proposal point X */
      U = (1-U)/(1-q);         /* U ~ U(0,1) */
      X = table[j] + U * d;    /* X ~ U(a_j, a_{j+1}) */

      if ( (j>=n) &&
	   ( (d>0 && X > table[n+3]) || (d<0 && X < table[n+3]) ) )
	/* reject points in terminal interval that are outside of slope */
	continue;

      /* acceptance/rejection step */
      V = _unur_call_urng(gen->urng);   /* V ~ U(0,1) */
      if ( V * (fmax-fmin) <= (PDF(X) - fmin) )  /* accept */
	return X;
      /* else try again */
    }

  } while(1);


} /* end of _unur_tabl_fast_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_tabl_fast_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /* (use fast table)                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
/*   UNUR_URNG *urng;             /\* pointer to uniform random number generator *\/ */
/*   double U,X,fx,V; */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TABL_fast_GEN,INFINITY);

  return 0.;

} /* end of _unur_tabl_rh_sample_check() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

#define empty_line() fprintf(log,"%s:\n",gen->genid);

void
_unur_tabl_fast_debug_init_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  double *table;
  int i,j,n;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_fast_GEN,RETURN_VOID);

  log = unur_get_stream();

  /* print slopes */
  fprintf(log,"%s: slopes for generator: %d\t(intervals: %d, area = %g)\n",gen->genid,
	  fastGEN->n_slopes, (int) fastGEN->n_ivs, fastGEN->bararea);
  table = fastGEN->table;
  n = 0;
  for (i=0; i<fastGEN->n_slopes; i++) {
    fprintf(log,"%s:  (%c) < %g, %g >\n",gen->genid,
	    (fastGEN->slopes[2*i]<fastGEN->slopes[2*i+1])?'-':'+',
	    fastGEN->slopes[2*i],fastGEN->slopes[2*i+1]);
    fprintf(log,"%s:\t[%2d] %g",gen->genid, (int)table[n], table[n+1]);
    for (j=n+2; j<=(int)table[n]; j++)
      fprintf(log,", %g", table[j]);
    fprintf(log,"  ||  %g, %g, %g\n",table[j],table[j+1],table[j+2]);
    n = j+3;
  }

  empty_line();
  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  empty_line();

} /* end of _unur_tabl_fast_debug_init_finished() */

/*****************************************************************************/

void
_unur_tabl_fast_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_fast_GEN,RETURN_VOID);

  log = unur_get_stream();

  empty_line();
  if (gen->status == UNUR_SUCCESS)
    fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  else
    fprintf(log,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
  empty_line();

  fflush(log);

} /* end of _unur_tabl_fast_debug_free() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
