/* ------------------------------------------------------------- */
/* File: example_vmt_1.c                                         */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

/* ------------------------------------------------------------- */

/* Example for method VMT with different marginal distributions. */
/*                                                               */
/* NOTICE that due to the transformation using Cholesky factors  */
/* of the covariance matrix the marginal distributions of the    */
/* resulting distribution are distorted.                         */
/* (For more details see manual for method VMT.)                 */

/* ------------------------------------------------------------- */

#define dim (3)

int main()
{
  int    i;
  double X[dim];

  /* Mean and covariance matrix for multivariate distribution.   */
  double mean[]  = { 1., 2., 3. };

  double covar[] = { 2., 2., 1., 
		     2., 4., 3.,
		     1., 3., 3. };

  /* parameters for beta "marginal" distribution.                */
  double beta_params[2];

  /* Declare the three UNURAN objects.                           */
  UNUR_DISTR *distr;    /* distribution object                   */
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  /* Create a multivariate distribution object with given mean,  */
  /* covariance matrix,and the "marginal" distributions          */
  /*   normal, Cauchy, and beta(3,5) distribution.               */
  distr = unur_distr_cvec_new(dim);
  unur_distr_cvec_set_mean(distr,mean);
  unur_distr_cvec_set_covar(distr,covar);
  beta_params[0] = 3; beta_params[1] = 5;
  unur_distr_cvec_set_stdmarginal_list(distr,
				       unur_distr_normal(NULL,0),
				       unur_distr_cauchy(NULL,0),
				       unur_distr_beta(beta_params,2));

  /* Choose a method: VMT.                                       */
  par = unur_vmt_new(distr);

  /* Create the generator object.                                */
  gen = unur_init(par);

  /* It is important to check if the creation of the generator   */
  /* object was successful. Otherwise `gen' is the NULL pointer  */ 
  /* and would cause a segmentation fault if used for sampling.  */
  if (gen == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }

  /* It is possible to reuse the distribution object to create   */
  /* another generator object. If you do not need it any more,   */
  /* it should be destroyed to free memory.                      */
  unur_distr_free(distr);

  /* Now you can use the generator object `gen' to sample from   */
  /* the distribution. Eg.:                                      */
  for (i=0; i<10; i++) {
    unur_sample_vec(gen, X);
    printf("(%f,%f,%f)\n", X[0], X[1], X[2]);
  }
 
  /* When you do not need the generator object any more, you     */
  /* can destroy it.                                             */
  unur_free(gen);

  exit (EXIT_SUCCESS);

} /* end of main() */

/* ------------------------------------------------------------- */
