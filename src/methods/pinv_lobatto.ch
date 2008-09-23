/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_lobatto.c                                               *
 *                                                                           *
 *   Routines for Gauss-Lobatto integration.                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/** Gauss-Lobatto integration                                               **/
/*****************************************************************************/

/* points for Lobatto integration */
#define W1 (0.17267316464601146)   /* = 0.5-sqrt(3/28) */
#define W2 (1.-W1)

/*---------------------------------------------------------------------------*/

double
_unur_pinv_lobatto5 (struct unur_gen *gen, double x, double h)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of the PDF over the interval (x,x+h)           */
     /* using Gauss-Lobatto integration with 5 points. (non-adaptive)        */
     /*                                                                      */
     /* Halfs intervals recursively if error is too large.                   */
     /*                                                                      */
     /* The recursion stops when the ABSOLUTE error is less than the         */
     /* respective given tolerances.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... left boundary point of interval                            */
     /*   h   ... length of interval                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{ 
  return (9*(PDF(x)+PDF(x+h))+49.*(PDF(x+h*W1)+PDF(x+h*W2))+64*PDF(x+h/2.))*h/180.;
} /* end of _unur_pinv_lobatto5() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_adaptivelobatto5 (struct unur_gen *gen, double x, double h, double tol,
			     struct unur_pinv_CDFtable *CDFtable)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of the PDF over the interval (x,x+h)           */
     /* using adaptive Gauss-Lobatto integration with 5 points.              */
     /*                                                                      */
     /* Halfs intervals recursively if error is too large.                   */
     /*                                                                      */
     /* The recursion stops when the ABSOLUTE error is less than the         */
     /* respective given tolerances.                                         */
     /*                                                                      */
     /* As a side effect, it stores boundaries and integrals for all         */
     /* subintervals in each recursion step where no further adaptation is   */
     /* required.                                                            */
     /* The values are stored in GEN->CDFtable (unless it equals NULL).      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   x        ... left boundary point of interval                       */
     /*   h        ... length of interval                                    */
     /*   tol      ... tolerated ABSOLUTE error                              */
     /*   CDFtable ... table for storing CDF values (may be NULL)            */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{ 
  double fl, fc, fr;  /* values of PDF at x, x+h/2, and x+h */
  double int1;        /* estimated values for integral */

  /* check length of interval */
  if (_unur_iszero(h))
    return 0.;

  /* arguments which are not finite (inf or NaN) cause infinite recursions */
  if (!_unur_isfinite(x+h)) {
    _unur_error(gen->genid,UNUR_ERR_INF,"boundaries of integration domain not finite");
    return INFINITY;
  }

  /* compute PDF values */
  fl = PDF(x);
  fc = PDF(x+h/2.);
  fr = PDF(x+h);
 
  /* first estimate for integral on [x,x+h] */
  int1 = (9*(fl+fr)+49.*(PDF(x+h*W1)+PDF(x+h*W2))+64*fc)*h/180.;

  /* run adaptive steps */
  return _unur_pinv_adaptivelobatto5_rec(gen,x,h,tol,int1,fl,fc,fr,CDFtable);

} /* end of _unur_pinv_adaptivelobatto5() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_adaptivelobatto5_rec (struct unur_gen *gen, double x, double h, double tol,
				 double int1, double fl, double fc, double fr,
				 struct unur_pinv_CDFtable *CDFtable)
     /*----------------------------------------------------------------------*/
     /* run recursion for adaptive Lobatto integration.                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   x        ... left boundary point of interval                       */
     /*   h        ... length of interval                                    */
     /*   tol      ... tolerated ABSOLUTE error                              */
     /*   fl       ... PDF at x                                              */
     /*   fc       ... PDF at x+h/2                                          */
     /*   fr       ... PDF at x+h                                            */
     /*   CDFtable ... table for storing CDF values (may be NULL)            */
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{
  double flc, frc;    /* values at PDF at x+h/4 and x+3*h/4 */
  double int2;        /* estimated values of integrals */
  double intl, intr;  /* left and right part of int2 */

  /* compute PDF values */
  flc = PDF(x+h/4);
  frc = PDF(x+3*h/4);
 
  /* compute integral on [x,x+h/2] and on [x+h/2,x+h] */
  intl = (9*(fl+fc)+49.*(PDF(x+h*W1*0.5)+PDF(x+h*W2*0.5))+64*flc)*h/360.;
  intr = (9*(fc+fr)+49.*(PDF(x+h*(0.5+W1*0.5))+PDF(x+h*(0.5+W2*0.5)))+64*frc)*h/360.;
  int2 = intl + intr;

  /* check whether accuracy goal is reached */
  if (fabs(int1-int2) < tol) 
    /* goal reached; nothing left to do */
    ;

  else {
    /* error above tolerance */
    if (_unur_FP_same(x+h/2.,x)) {
      /* we cannot decrease length of subintervals any more */
      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
		    "numeric integration did not reach full accuracy");
      /* Remark: Since we are halving intervals, this comparision */
      /* limits the maximal number of iterations to at most 2048. */
    }
    else {
      /* recompute with shorter intervals */
      return ( _unur_pinv_adaptivelobatto5_rec(gen,x,    h/2,tol,intl,fl,flc,fc,CDFtable) +
	       _unur_pinv_adaptivelobatto5_rec(gen,x+h/2,h/2,tol,intr,fc,frc,fr,CDFtable) );
    }
  }

#ifdef PINV_USE_CDFTABLE
  /* store integral values */
  if (CDFtable) {
    /* l.h.s. subinterval */
    _unur_pinv_CDFtable_append(CDFtable, x+h/2., intl);
    /* r.h.s. subinterval */
    _unur_pinv_CDFtable_append(CDFtable, x+h, intr);
  }
  /* Remark: we do not throw a warning if the table size is exceeded. */
#endif

  /* return estimate for integral */
  return int2;

} /* end of _unur_pinv_adaptivelobatto5_rec() */

/*---------------------------------------------------------------------------*/
#undef W1
#undef W2
/*---------------------------------------------------------------------------*/

double
_unur_pinv_Udiff_lobatto (struct unur_gen *gen, double x, double h, double utol)
     /*----------------------------------------------------------------------*/
     /* Compute difference CDF(x+h)-CDF(x) (approximately), where CDF is the */
     /* integral of the given (quasi-) density.                              */
     /* It makes use the the table of CDF values computed in                 */
     /* _unur_pinv_pdfarea().                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... left boundary point of interval                            */
     /*   h   ... length of interval                                         */
     /*   tol ... tolerated ABSOLUTE error                                   */
     /*                                                                      */
     /* return:                                                              */
     /*    (approximate) difference CDF(x+h) - CDF(x)                        */
     /*----------------------------------------------------------------------*/
{
#ifdef PINV_USE_CDFTABLE
  struct unur_pinv_CDFtable *CDFtable = GEN->CDFtable; /* table of CDF values */
  int cur;                    /* pointer to current position in table */
  double x1;                  /* left and right boundary of a subinterval */
  double Udiff;               /* value to be returned */
#endif

  /* arguments which are not finite (inf or NaN) cause infinite recursions */
  if (!_unur_isfinite(x+h)) {
    /* _unur_warning(gen->genid,UNUR_ERR_INF,"boundaries of integration domain not finite"); */
    return INFINITY;
  }


#ifdef PINV_USE_CDFTABLE
  /* check for table */
  if (CDFtable == NULL) {
    /* there is no table: use adaptive integration */
    
    if (x < GEN->bleft || x+h > GEN->bright)
      /* the interval [x,x+h] is too large. thus we use simple Lobatto integration. */
      return _unur_pinv_lobatto5(gen, x, h);
    else
      return _unur_pinv_adaptivelobatto5(gen, x, h, utol, NULL);
  }

  /* else: we try to read CDF values from table */

  /* move pointer for reading to start position in interval */
  cur = CDFtable->cur_iv;

  /* first entry in interval */
  while (cur < CDFtable->n_values &&
	 CDFtable->values[cur].x < x)
    ++cur;

  /* did we find such an entry ? */
  if (cur >= CDFtable->n_values) {
    /* we must use adaptive Lobatto integration if the table for  */
    /* CDF values was too small.                                  */
    return _unur_pinv_adaptivelobatto5(gen, x, h, utol, NULL);
  }

  /* store x value and goto next entry */
  x1 = CDFtable->values[cur].x;
  ++cur;

  /* are there more than one entry in interval ? */
  if (cur >= CDFtable->n_values ||
      CDFtable->values[cur].x > x+h) {
    /* there is at most one entry in the interval [x,x+h]. */
    /* thus we (can) use simple Lobatto integration.       */
    return _unur_pinv_lobatto5(gen, x, h);
  }

  /* there are more than one entries in the interval.              */
  /* thus there is at least one subinterval from adapative Lobatto */
  /* integration within [x,x+h]. we reuse the stored values.       */
  /* for the remaining two subintervals at the boundary [x,x+h]    */
  /* we use simple Lobatto integration.                            */

  Udiff = _unur_pinv_lobatto5(gen, x, x1 - x);
  do {
    Udiff += CDFtable->values[cur].u;
    /*       = _unur_pinv_lobatto5(gen, x1, x2 - x1) */
    x1 = CDFtable->values[cur].x;
    ++cur;
  } while (cur < CDFtable->n_values
	   && CDFtable->values[cur].x <= x+h);

  /* We have to distinguish two cases: */
  if (x+h < GEN->bright && cur >= CDFtable->n_values) {
    /* the table of CDF values is too small */
    Udiff += _unur_pinv_adaptivelobatto5(gen, x1, x+h, utol, NULL);
  }
  else {
    /* the table is not too small but x+h is outside the computational domain */
    Udiff += _unur_pinv_lobatto5(gen, x1, x+h - x1);
  }
      
  return Udiff;

#else

#ifdef PINV_USE_SIMPLE_LOBATTO
  return _unur_pinv_lobatto5(gen, x, h);
#else
  if (x < GEN->bleft || x+h > GEN->bright)
    /* the interval [x,x+h] is too large. thus we use simple Lobatto integration. */
    return _unur_pinv_lobatto5(gen, x, h);
  else
    return _unur_pinv_adaptivelobatto5(gen, x, h, utol, NULL);
#endif

#endif /* defined(PINV_USE_CDFTABLE) */

} /* end of _unur_pinv_Udiff() */


/*---------------------------------------------------------------------------*/
#ifdef PINV_USE_CDFTABLE
/*---------------------------------------------------------------------------*/

struct unur_pinv_CDFtable *
_unur_pinv_CDFtable_create (int size)
     /*----------------------------------------------------------------------*/
     /* create table of CDF values.                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   size ... size of table                                             */
     /*----------------------------------------------------------------------*/
{
  struct unur_pinv_CDFtable *table;

  /* check argument */
  if (size<=0)
    return NULL;

  /* allocate memory */
  table = _unur_xmalloc( sizeof(struct unur_pinv_CDFtable) );
  table->values = _unur_xmalloc(size * sizeof(struct unur_pinv_CDFvalues) ); 

  /* set counter */
  table->size = size;
  table->n_values = 0;
  table->cur_iv = 0;

  return table;
} /* end of _unur_pinv_CDFtable_create() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_CDFtable_append (struct unur_pinv_CDFtable *table, double x, double u)
     /*----------------------------------------------------------------------*/
     /* append entry to table of CDF values.                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   tables ... table with CDF values                                   */
     /*   x      ... right boundary of subinterval                           */
     /*   u      ... integral over subinterval                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  if (table==NULL) 
    return UNUR_ERR_NULL;

  if (table->n_values >= table->size - 1)
    /* we do not write a warning here */
    return UNUR_ERR_GENERIC;

  table->values[table->n_values].x = x;
  table->values[table->n_values].u = u;
  ++(table->n_values);

  return UNUR_SUCCESS;
} /* end of _unur_pinv_CDFtable_append() */

/*---------------------------------------------------------------------------*/

void 
_unur_pinv_CDFtable_resize (struct unur_pinv_CDFtable **table)
     /*----------------------------------------------------------------------*/
     /* resize table of CDF values.                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   table ... pointer to pointer to table with CDF values              */
     /*----------------------------------------------------------------------*/
{
  if (*table) {
    *table = _unur_xrealloc(*table, (*table)->n_values * sizeof(double));
    (*table)->size = (*table)->n_values;
  }
  /* else: nothing to do */
  
} /* end of _unur_pinv_CDFtable_resize() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_CDFtable_free (struct unur_pinv_CDFtable **table)
     /*----------------------------------------------------------------------*/
     /* destroy table of CDF values and set pointer to NULL.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   table ... pointer to pointer to table with CDF values              */
     /*----------------------------------------------------------------------*/
{
  if (*table) {
    free ((*table)->values);
    free (*table);
    *table = NULL;
  }
  /* else: nothing to do */

} /* end of _unur_pinv_CDFtable_free() */

/*---------------------------------------------------------------------------*/
#endif /* defined(PINV_USE_CDFTABLE) */
/*---------------------------------------------------------------------------*/
