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

#include <source_unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

/* defaults */
#define CHI2_CLASSMIN_DEFAULT  20  /* minimum number of observations in class */
#define CHI2_INTERVALS_DEFAULT 50  /* number of intervals for chi^2 test     */

/* constant */
#define CHI2_SIZE_PMF_VECTOR   100 /* size of probability vector when 
					only p.m.f. is given.  */
/** TODO: make more flexible **/

/*---------------------------------------------------------------------------*/
static char test_name[] = "Chi^2-Test";
/*---------------------------------------------------------------------------*/

static double _unur_test_chi2_discr( struct unur_gen *gen, int samplesize, int classmin, int verbose );
static double _unur_test_chi2_cont( struct unur_gen *gen, int intervals, int samplesize, int classmin, int verbose );

/*---------------------------------------------------------------------------*/

double
unur_test_chi2( struct unur_gen *gen, 
		int intervals,
		int samplesize, 
		int classmin,
		int verbose )
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
     /*   verbose  ... verbosity level                                       */
     /*                  0 = no output on stdout                             */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1.                                                         */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL(test_name,gen,-1.);

  if (verbose >= 1)
    printf("\nGOODNESS-OF-FIT TESTS:\n");

  switch (gen->method & UNUR_MASK_TYPE) {

  case UNUR_METH_DISCR:
    return _unur_test_chi2_discr(gen, samplesize, classmin, verbose);

  case UNUR_METH_CONT:
    return _unur_test_chi2_cont(gen, intervals, samplesize, classmin, verbose);

  case UNUR_METH_VEC:
    _unur_error(test_name,UNUR_ERR_GENERIC,"Not implemented for multivariate distributions!");
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
		       int verbose )
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
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1.                                                         */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr.data.discr
  double *prob;         /* pointer to probability vectors */
  int n_prob;           /* length of probability vector   */

  int *observed;        /* vector for observed occurrences */
  double pval;          /* p-value */
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,-1.);

  /* check distribution: discrete univariate */
  if ( (gen->method & UNUR_MASK_TYPE) != UNUR_METH_DISCR ) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"wrong distribution type");
    return -1.;
  }

  /* pointer to distribution object */
  prob = DISTR.prob;
  n_prob = DISTR.n_prob;

#if 0
  /* check argument: need probability vector */
  if (prob == NULL) {
    /* no probability vector  -->  p.m.f. required */
    if (DISTR.pmf == NULL) {
      _unur_error(test_name,UNUR_ERR_GENERIC,"probability vector or p.m.f. required");
      return -1.;
    }
    /* have to make own probability vector */
    n_prob = CHI2_SIZE_PMF_VECTOR;
    prob = _unur_malloc( n_prob * sizeof(double) );
    if (!prob) return -1.;
    for (i=0; i<n_prob; i++) {
      /** TODO: here we assume domain [0,infinity] */
      /** furthermore we cut domain at CHI2_SIZE_PMF_VECTOR */
      prob[i] = _unur_discr_PMF( i, &(gen->distr) );
    }
  }
#endif
  /* check argument: need probability vector */
  if (prob == NULL)
    return -1.;

  /* allocate memory for observations */
  observed = _unur_malloc( n_prob * sizeof(int));

  /* clear array */
  for( i=0; i<n_prob; i++ )
    observed[i] = 0;

  /* samplesize */
  if( samplesize <= 0 )
    samplesize = (INT_MAX/n_prob > n_prob) ? n_prob*n_prob : INT_MAX;

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    /* sample */
    j = _unur_sample_discr(gen);
    /* shift vector */
    j -= DISTR.domain[0];
    /* check range of random variates !! */
    if (j >= n_prob)
      j = n_prob - 1;  /* put into the last bin */
    ++observed[j];
  }

  if (verbose >= 1) {
    printf("\nChi^2-Test for discrete distribution with given probability vector:");
    printf("\n  length     = %d\n",n_prob);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(prob, observed, n_prob, classmin, verbose);

  /* free memory */
  free(observed);
  if (DISTR.prob == NULL) free(prob);

  /* return result of test */
  return pval;

#undef DISTR
} /* end of _unur_test_chi2_discr() */

/*---------------------------------------------------------------------------*/

static double
_unur_test_chi2_cont(struct unur_gen *gen, 
		     int intervals, 
		     int samplesize, 
		     int classmin,
		     int verbose )
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
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* comment:                                                             */
     /*   uses routine chi2test().                                           */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr.data.cont

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
    return -1.;
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

  /* compute Fl and Fr */
  if (gen->distr.set & UNUR_DISTR_SET_TRUNCATED) {
    Fl = (DISTR.trunc[0] <= -INFINITY) ? 0. : cdf(DISTR.trunc[0],&(gen->distr));
    Fr = (DISTR.trunc[1] >=  INFINITY) ? 1. : cdf(DISTR.trunc[1],&(gen->distr));
  }
  else {
    Fl = (DISTR.domain[0] <= -INFINITY) ? 0. : cdf(DISTR.domain[0],&(gen->distr));
    Fr = (DISTR.domain[1] >=  INFINITY) ? 1. : cdf(DISTR.domain[1],&(gen->distr));
  }
  Fdelta = Fr - Fl;

  /* Fr - Fl <= 0. is a fatal error */
  if (Fdelta <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
    return -1.;
  }

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    F = cdf( _unur_sample_cont(gen), &(gen->distr) );
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
    printf("\nChi^2-Test for continuous distribution:");
    printf("\n  intervals  = %d\n",intervals);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(NULL, observed, intervals, classmin, verbose );

  /* free memory */
  free(observed);

  /* return result of test */
  return pval;

#undef DISTR
} /* end of _unur_test_chi2_cont() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2test( double *prob, 
		     int *observed, 
		     int len, int classmin, int verbose )
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
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1.                                                         */
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

  static UNUR_DISTR *chisquare_distr = NULL; /* distribution object to chi^2 distribution */

  /* check arguments */
  CHECK_NULL(observed,-1.);

  /* make distribution object for chi^2 distribution */
  if (chisquare_distr == NULL) {
    df = 1.;  /* dummy value for degrees of freedom (will be set later) */
    chisquare_distr = unur_distr_chisquare( &df, 1 );
  }

  /* do we have a CDF for the chi^2 distribution ? */
  if (chisquare_distr->data.cont.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF for CHI^2 distribution required");
    return -1.;
  }

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
	printf("Class #%d:\tobserved %d\texpected %.2f\n",classes,clobsd,clexpd);
      clexpd = 0.;
      clobsd = 0;
      classes++;
    }
  }

  /* there must be at least two classes */
  if (classes < 2) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"too few classes!");
    if (verbose >= 1)
      printf("\nCannot run chi^2-Test: too few classes\n");
    return -1.;
  }

  /* evaluate test statistics */
  /* set degrees of freedom */
  df = (double)(classes-1);
  unur_distr_cont_set_pdfparams( chisquare_distr, &df, 1 );
  /* p-value */
  pval = 1. - _unur_cont_CDF( chi2, chisquare_distr );

  /* print result (if requested) */
  if (verbose >= 1) {
    printf("\nResult of chi^2-Test:\n  samplesize = %d\n",samplesize);
    printf("  classes    = %d\t (minimum per class = %d)\n", classes, classmin);
    printf("  chi2-value = %g\n  p-value    = %g\n\n", chi2, pval);
  }

  /* return result of test */
  return pval;

} /* end of _unur_test_chi2test() */

/*---------------------------------------------------------------------------*/
