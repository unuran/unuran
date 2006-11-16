/* ------------------------------------------------------------- */
/* File: example_sample_wo_replacement.c                         */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

/* ------------------------------------------------------------- */

int sample_wo_replacement( int *s, int k, double *pv, int n );
/* draw a sample of size k from {0,...,n-1} without replacement    */
/* with probability vector pv                                    */

/* ------------------------------------------------------------- */

int sample_wo_replacement( int *s, int k, double *pv, int n )
{
  /* Declare the three UNURAN objects.                           */
  UNUR_DISTR *distr;    /* distribution object                   */
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  double sum, *gpv;
  int i;

  /* check arguments */
  if (k<1 || k>n || n < 1) return UNUR_FAILURE;

  /* create discrete distribution object */
  distr = unur_distr_discr_new();
  if (unur_distr_discr_set_pv(distr, pv, n) != UNUR_SUCCESS)
    return UNUR_FAILURE;
  for (sum=0., i=0; i<n; i++) sum += pv[i];
  unur_distr_discr_set_pmfsum(distr,sum);

  /* create generator object */
  par = unur_dss_new(distr);
  unur_set_use_distr_privatecopy(par,FALSE);
  gen = unur_init(par);

  /* we need the pointer to the internal PV    */
  /* Remark: this produces a compiler warning: */
  /*   'incompatible pointer type'             */
  unur_distr_discr_get_pv(distr,&gpv);

  /* draw a sample */
  for (i=0; i<k; i++) {
    s[i] = unur_sample_discr(gen);
    sum -= gpv[s[i]];
    unur_distr_discr_set_pmfsum(distr,sum);
    gpv[s[i]] = 0.;
    /* unur_reinit(gen); ... not required for this simple method */
  }

  /* free memory */
  unur_distr_free(distr);
  unur_free(gen);

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of sample_wo_replacement() */

/* ------------------------------------------------------------- */

int main(void)
{
#define N  (10)
  int i,j;
  double pv[N] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
  int n = N;

  int s[N];
  int k = 8;

  for (j=0; j<10; j++) {
    sample_wo_replacement( s, k, pv, n );
    for (i=0; i<k; i++) printf("%d, ",s[i]);
    printf("\n");
  }

  exit (EXIT_SUCCESS);

} /* end of main() */

/* ------------------------------------------------------------- */

