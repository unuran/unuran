/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      chi2test.c                                                   *
 *                                                                           *
 *   Chi^2 Tests (to validate implementation)                                *
 *                                                                           *
 *   WARNING!                                                                *
 *   A succesfull chi^2 test shows that the implementation is (probably)     *
 *   correct. It DOES NOT mean that the generator is of good quality!        *
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

#include <limits.h>
#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include <distr/discr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <specfunct/unur_specfunct_source.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/

/* defaults */
#define CHI2_CLASSMIN_DEFAULT  20 /* minimum number of observations in class */
#define CHI2_SAMPLEFAC  40  
         /* if samplesize<=0 use samplesize = CHI2_SAMPLEFAC * intervals^dim */
#define CHI2_INTERVALS_DEFAULT 50  /* number of intervals for chi^2 test     */
#define CHI2_MAX_DIMENSIONS 40 /* max number of dimensions for CHI2 tests  */

/* constants */
#define CHI2_MAX_SAMPLESIZE 1000000

/*---------------------------------------------------------------------------*/
static char test_name[] = "Chi^2-Test";
/*---------------------------------------------------------------------------*/

static double _unur_test_chi2_discr( struct unur_gen *gen, int samplesize, int classmin,
				     int verbose, FILE *out );

static double _unur_test_chi2_cont( struct unur_gen *gen, int intervals, int samplesize, int classmin,
				    int verbose, FILE *out );

static double _unur_test_chi2_cemp( struct unur_gen *gen, int intervals, int samplesize, int classmin,
				    int verbose, FILE *out );

static double _unur_test_chi2_vec( struct unur_gen *gen, int intervals, int samplesize, int classmin,
                                    int verbose, FILE *out );

static double _unur_test_chi2test( double *prob, int *observed, int len, int classmin,
				   int verbose, FILE *out );

/*---------------------------------------------------------------------------*/

double
unur_test_chi2( struct unur_gen *gen, 
		int intervals,
		int samplesize, 
		int classmin,
		int verbose,
		FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for univariate distributions                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator                                */
     /*   intervals  ... number if intervals                                 */
     /*                  case probability vector: use its length             */
     /*                  otherwise and if <= 0 use default                   */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, (#classes)^2 is used as default.)         */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output on out                                */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL(test_name,gen,-1.);

  if (verbose >= 1)
    fprintf(out,"\nGOODNESS-OF-FIT TESTS:\n");

  switch (gen->method & UNUR_MASK_TYPE) {

  case UNUR_METH_DISCR:
    return _unur_test_chi2_discr(gen, samplesize, classmin, verbose, out);

  case UNUR_METH_CONT:
    return _unur_test_chi2_cont(gen, intervals, samplesize, classmin, verbose, out);

  case UNUR_METH_CEMP:
    return _unur_test_chi2_cemp(gen, intervals, samplesize, classmin, verbose, out);

  case UNUR_METH_VEC:
    /* TODO : testing for other than normal distributions */
    return _unur_test_chi2_vec(gen, intervals, samplesize, classmin, verbose, out);

  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"Not implemented for such distributions!");
    return -1.;
  }

  /* this statement should be reached */
  return -1.;
} /* end of unur_test_chi2() */

/*---------------------------------------------------------------------------*/

static double
_unur_test_chi2_discr( struct unur_gen *gen, 
		       int samplesize, 
		       int classmin, 
		       int verbose,
		       FILE *out)
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for univariate discrete distributions                     */
     /* with given probability vector.                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, (#classes)^2 is used as default.)         */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output on stdout                             */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.discr
  double *pv;           /* pointer to probability vectors */
  int n_pv;             /* length of probability vector   */

  int *observed;        /* vector for observed occurrences */
  double pval;          /* p-value */
  int had_PV;           /* whether we had a PV for test or not */
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,-1.);

  /* check distribution: discrete univariate */
  if ( (gen->method & UNUR_MASK_TYPE) != UNUR_METH_DISCR ) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"wrong distribution type");
    return -1.;
  }

  /* probability vector */
  if (DISTR.pv == NULL) {
    had_PV = FALSE;
    /* no PV given --> try to compute PV */
    if (!unur_distr_discr_make_pv( gen->distr )) {
      /* not successful */
      return -2.;
    }
  }
  else 
    had_PV = TRUE;
  /* pointer to PV */
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;

  /* allocate memory for observations */
  observed = _unur_malloc( n_pv * sizeof(int));

  /* clear array */
  for( i=0; i<n_pv; i++ )
    observed[i] = 0;

  /* samplesize */
  if( samplesize <= 0 ) {
    samplesize = (INT_MAX/n_pv > n_pv) ? n_pv * n_pv : 1000000;
    samplesize = max(samplesize,1000000);
  }

  samplesize = min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    /* sample */
    j = _unur_sample_discr(gen);
    if (verbose >= 3) fprintf(out,"i = %d\n",j);
    /* shift vector */
    j -= DISTR.domain[0];
    /* check range of random variates !! */
    if (j < n_pv)
      ++observed[j];
    /* else: ignore number --> chop off tail */
  }

  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for discrete distribution with given probability vector:");
    fprintf(out,"\n  length     = %d\n",n_pv);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(pv, observed, n_pv, classmin, verbose, out);

  /* free memory */
  free(observed);
  if (!had_PV) {
    free (DISTR.pv);
    DISTR.pv = NULL;
    DISTR.n_pv = 0;
  }

  /* return result of test */
  return pval;

#undef DISTR
} /* end of _unur_test_chi2_discr() */

/*---------------------------------------------------------------------------*/

static double
_unur_test_chi2_cont( struct unur_gen *gen, 
		      int intervals, 
		      int samplesize, 
		      int classmin,
		      int verbose,
		      FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for univariate continuous distributions.                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   intervals  ... number of intervals in which (0,1) is partitioned   */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, intervals^2 is used as default.)          */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output on stdout                             */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.cont

  double F, Fl, Fr, Fdelta;  /* value of CDF (at left and right boundary point) */
  UNUR_FUNCT_CONT *cdf;      /* pointer to CDF */
  int *observed;             /* vector for observed occurrences */
  double pval;               /* p-value */
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,-1.);
  /* we do not check magic cookies here */

  /* CDF required */
  cdf = DISTR.cdf;
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required for continuous random variates!");
    return -2.;
  }

  /* check given number of intervals */
  if (intervals <= 2)
    intervals = CHI2_INTERVALS_DEFAULT;

  /* allocate memory for observations */
  observed = _unur_malloc( intervals * sizeof(int));

  /* clear array */
  for( i=0; i<intervals; i++ )
    observed[i] = 0;

  /* samplesize */
  if( samplesize <= 0 )
    samplesize = (INT_MAX/intervals > intervals) ? intervals*intervals : INT_MAX;

  samplesize = min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* compute Fl and Fr */
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    Fl = (DISTR.trunc[0] <= -INFINITY) ? 0. : cdf(DISTR.trunc[0], gen->distr);
    Fr = (DISTR.trunc[1] >=  INFINITY) ? 1. : cdf(DISTR.trunc[1], gen->distr);
  }
  else {
    Fl = (DISTR.domain[0] <= -INFINITY) ? 0. : cdf(DISTR.domain[0], gen->distr);
    Fr = (DISTR.domain[1] >=  INFINITY) ? 1. : cdf(DISTR.domain[1], gen->distr);
  }
  Fdelta = Fr - Fl;

  /* Fr - Fl <= 0. is a fatal error */
  if (Fdelta <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
    free (observed);
    return -1.;
  }

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    if (verbose >= 3) {
      double x = _unur_sample_cont(gen);
      F = cdf( x, gen->distr );
      fprintf(out,"x = %g\n",x);
    }
    else {
      F = cdf( _unur_sample_cont(gen), gen->distr );
    }
    F = (F-Fl)/Fdelta;
    j = (int)(intervals * F);
    if (j > intervals) {   
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) > Fmax (out of domain).");
      j = intervals-1;
    }
    if (j >= intervals)    /* cdf() might return 1. */
      j = intervals-1;
    if (j < 0 ) {           /* there is something wrong with the boundaries */
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) < 0 (out of domain).");
      j = 0;
    }
    ++observed[j];
  }

  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for continuous distribution:");
    fprintf(out,"\n  intervals  = %d\n",intervals);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(NULL, observed, intervals, classmin, verbose, out );

  /* free memory */
  free(observed);

  /* return result of test */
  return pval;

#undef DISTR
} /* end of _unur_test_chi2_cont() */

/*---------------------------------------------------------------------------*/

static double
_unur_test_chi2_cemp( struct unur_gen *gen, 
		      int intervals, 
		      int samplesize, 
		      int classmin,
		      int verbose,
		      FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for continuous empirical distributions.                   */
     /* Tests for standard normal distribution only!                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   intervals  ... number of intervals in which (0,1) is partitioned   */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, intervals^2 is used as default.)          */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output on stdout                             */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
  UNUR_DISTR *distr_normal;  /* object for standard normal distribution */
  UNUR_FUNCT_CONT *cdf;      /* pointer to CDF */
  double F;                  /* value of CDF */
  int *observed;             /* vector for observed occurrences */
  double pval;               /* p-value */
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,-1.);
  /* we do not check magic cookies here */

  /* CDF for standard normal distribution */
  distr_normal = unur_distr_normal( NULL, 0 );
  cdf = distr_normal->data.cont.cdf;

  /* check given number of intervals */
  if (intervals <= 2)
    intervals = CHI2_INTERVALS_DEFAULT;

  /* allocate memory for observations */
  observed = _unur_malloc( intervals * sizeof(int));

  /* clear array */
  for( i=0; i<intervals; i++ )
    observed[i] = 0;

  /* samplesize */
  if( samplesize <= 0 )
    samplesize = (INT_MAX/intervals > intervals) ? intervals*intervals : INT_MAX;
  samplesize = min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    F = cdf( _unur_sample_cont(gen), distr_normal );
    j = (int)(intervals * F);
    if (j > intervals) {   
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) > Fmax (out of domain).");
      j = intervals-1;
    }
    if (j >= intervals)    /* cdf() might return 1. */
      j = intervals-1;
    if (j < 0 ) {           /* there is something wrong with the boundaries */
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) < 0 (out of domain).");
      j = 0;
    }
    ++observed[j];
  }

  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for continuous empirical distribution:");
    fprintf(out,"\n(Assumes standard normal distribution!)");
    fprintf(out,"\n  intervals  = %d\n",intervals);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(NULL, observed, intervals, classmin, verbose, out );

  /* free memory */
  _unur_distr_free(distr_normal);
  free(observed);

  /* return result of test */
  return pval;

} /* end of _unur_test_chi2_cemp() */

/*---------------------------------------------------------------------------*/

static double
_unur_test_chi2_vec ( struct unur_gen *gen, 
		      int intervals, 
		      int samplesize, 
		      int classmin,
		      int verbose,
		      FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for multivariate (NORMAL) continuous distributions.       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   intervals  ... number of intervals in which (0,1) is partitioned   */
     /*                  (each dimension is partitioned in #intervals)       */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, CHI2_SAMPLEFAC * intervals^dim is used)   */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output on stdout                             */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{

  int dim;         /* dimension of multivariate distribution */
  double *z;       /* sampling vector */
  int *bg[CHI2_MAX_DIMENSIONS]; /* vectors for observed occurrences */
  double pval;     /* p-value */
  int i,j, sumintervals;
  int *idx;	   /* index array */
  int dimintervals[CHI2_MAX_DIMENSIONS]; /* for marginal chi2 tests */ 
  int totalintervals; /* sum of the dimintervals[] */
  
  /* check arguments */
  CHECK_NULL(gen,-1.);
  /* we do not check magic cookies here */

  /* check given number of intervals */
  if (intervals <= 2)
    intervals = CHI2_INTERVALS_DEFAULT;

  dim = gen->distr->dim;
  if (dim < 2) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"Distribution dimension < 2 ?");
    return -1; 
  }

  if (dim > CHI2_MAX_DIMENSIONS) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"Distribution dimension too large");
    return -1; 
  }

  /* setup of intervals for each dimension */
  totalintervals=0;
  for (i=0; i<dim; i++) {
    dimintervals[i] = (int) ( intervals * (1.+UNUR_EPSILON)/(1+i) ) ;  
    if (dimintervals[i]>intervals) dimintervals[i] = intervals; /* in case we wish to make changes in previous line */
    if (dimintervals[i]<2) dimintervals[i]=2; /* 1 would be safer, but makes no sense in this context */
    totalintervals += dimintervals[i];
  }

  totalintervals += 0 ; /* maybe we'll need this parameter someday */

  /* allocate memory */
  idx = _unur_malloc( dim * sizeof(int));
  if (idx==NULL) {
      _unur_error(test_name,UNUR_ERR_GENERIC,"malloc error : idx");
      pval=-1; goto free_memory;
  }

  z = _unur_malloc( dim * sizeof(double));
  if (z==NULL) {
     _unur_error(test_name,UNUR_ERR_GENERIC,"malloc error : z");
     pval=-1; goto free_memory;
  }
  
  for (i=0; i<dim; i++) bg[i]=NULL;
  
  for (i=0; i<dim; i++) {
  bg[i] = _unur_malloc( dimintervals[i] * sizeof(int));
  if (bg[i]==NULL) {
     _unur_error(test_name,UNUR_ERR_GENERIC,"malloc error : bg");
     pval=-1; goto free_memory;
  }}   

  /* clear arrays */

  for (i=0; i<dim; i++) (void) memset(bg[i] , 0, dimintervals[i] * sizeof(int));
  (void) memset(idx, 0, dim * sizeof(int));

  /* samplesize */
  if( samplesize <= 0 ) {
    samplesize = abs(samplesize);
  }
  samplesize = min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    _unur_sample_vec(gen, z);
    sumintervals=0;
    for (j=0; j<dim; j++) {
      idx[j] = (int)( dimintervals[j] * _unur_sf_cdfnormal(z[j]) );

      if (idx[j]==dimintervals[j]) idx[j]--; /* cdf can return 1 ? */
      bg[j][idx[j]] += 1 ;
      sumintervals += dimintervals[j];
    }
  }

  /* and now make chi^2 test (marginal) */
  sumintervals = 0;
  for (j=0; j<dim; j++) {
    if (verbose >= 1) {
      fprintf(out,"\nChi^2-Test for multivariate continuous distribution\n");
      fprintf(out,"  samplesize = %d\n",samplesize);
      fprintf(out,"  intervals  = %d\n",dimintervals[j]);
      fprintf(out,"  marginal   = %d\n",j);

      for (i=0; i<dimintervals[j]; i++) {
        fprintf(out,"  bg[%d][%d]     = %d\n",j ,i ,bg[j][i] );
      }

    }

    pval = _unur_test_chi2test(NULL, bg[j] , dimintervals[j], classmin, verbose, out );
    sumintervals += dimintervals[j];
    
  }


free_memory:
  /* free memory */
  (idx==NULL) ? : free(idx);
  (z==NULL) ? : free(z);
  for (i=0; i<dim; i++) (bg[i]==NULL) ? : free(bg[i]);

  /* return result of test */

  return pval;

} /* end of _unur_test_chi2_vec() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2test( double *prob, 
		     int *observed, 
		     int len,
		     int classmin,
		     int verbose,
		     FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for discrete distributions.                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prob     ... probability vector (need not sum to 1)                */
     /*                  NULL indicates a uniform distribution on (0,len-1)  */
     /*   observed ... vector containing the observed number of occurrences  */
     /*   len      ... length of (both) vectors                              */
     /*   classmin ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose  ... verbosity level                                       */
     /*                  0 = no output on stdout                             */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out      ... output stream                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*                                                                      */
     /* comment:                                                             */
     /*   joines classes if expected number of occurrences in some classes   */
     /*   is too small, i.e., less than "classmin".                          */
     /*----------------------------------------------------------------------*/
{
  double chi2 = 0.;     /* chi2-value */
  double df;            /* degrees of freedom for chi^2 distribution */
  double pval;          /* p-value */
  double clexpd = 0.;   /* expected occurrences in class */
  int clobsd = 0;       /* observed occurrences in class */
  int classes = 0;      /* (total) number of classes     */
  double probsum = 0.;  /* sum of all "probabilities" (need not be 1, for convenience) */
  int samplesize = 0;
  double factor;        /* factor for calculating expected number of occurrences */
  int i;

  UNUR_DISTR *distr_chisquare = NULL; /* distribution object to chi^2 distribution */

  /* check arguments */
  CHECK_NULL(observed,-1.);

  /* minimum number of occurrences in a class */
  classmin = (classmin > 0) ? classmin : CHI2_CLASSMIN_DEFAULT;

  /* compute sample size */
  for( i=0; i<len; i++ )
    samplesize += observed[i];

  /* sum of probabilities (if not uniform distribution) */
  if (prob != NULL) {
    for( i=0; i<len; i++ )
      probsum += prob[i];
    factor = samplesize/probsum;
  }
  else   /* uniform distribution on (0,len-1) --> prob[i] = 1/len */
    factor = ((double)samplesize)/len;

  /* compute chi^2 value */
  for( i=0; i<len; i++ ) {
    clexpd += (prob) ? prob[i]*factor : factor;  /* expected occurrences in this class */
    clobsd += observed[i];                       /* observed occurrences in this class */
    
    if (clexpd >= classmin || i == len-1) {
      /* number of expected occurrences is large enough or end of array */
      chi2 += (clobsd-clexpd)*(clobsd-clexpd)/clexpd;
      if (verbose >= 2)
	fprintf(out,"Class #%d:\tobserved %d\texpected %.2f\n",classes,clobsd,clexpd);
      clexpd = 0.;
      clobsd = 0;
      classes++;
    }
  }

  /* there must be at least two classes */
  if (classes < 2) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"too few classes!");
    if (verbose >= 1)
      fprintf(out,"\nCannot run chi^2-Test: too few classes\n");
    return -1.;
  }

  /* evaluate test statistics */
  /* chisquare distribution with df degrees of freedom */
  df = (double)(classes-1);
  distr_chisquare = unur_distr_chisquare( &df, 1 );
  /* p-value */
  if (distr_chisquare->data.cont.cdf) {
    pval = 1. - _unur_cont_CDF( chi2, distr_chisquare );
  }
  else {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF for CHI^2 distribution required");
    pval = -2.;
  }
  _unur_distr_free(distr_chisquare);

  /* print result (if requested) */
  if (verbose >= 1 && pval >= 0.) {
    fprintf(out,"\nResult of chi^2-Test:\n  samplesize = %d\n",samplesize);
    fprintf(out,"  classes    = %d\t (minimum per class = %d)\n", classes, classmin);
    fprintf(out,"  chi2-value = %g\n  p-value    = %g\n\n", chi2, pval);
  }

  /* return result of test */
  return pval;

} /* end of _unur_test_chi2test() */

/*---------------------------------------------------------------------------*/
