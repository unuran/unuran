/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Common test routines                                                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "testunuran.h"

/*---------------------------------------------------------------------------*/

static struct prng *urng = NULL;

/*---------------------------------------------------------------------------*/

static struct list_distr *list_of_distr = NULL;  /* pointer to list of distributions */
static int n_distr;                       /* number of distributions in list */

#define N_DISTRIBUTIONS 1000         /* maximum number of distributions       */

/*---------------------------------------------------------------------------*/

static double *make_random_vector(int len);
/* make a probability vector (use random entries) */

static double *make_geometric_vector(double q, int len); 
/* make a probability vector (use geometric distribution) */

/*---------------------------------------------------------------------------*/

/* make list of distributions                                                */
int make_list_of_distributions( struct list_distr **return_list )
{
  double fpar[10];
  double *prob;
  struct list_distr *list;
  UNUR_DISTR *distr;

  /* we use Mersenne Twister as uniform random number generator */
  if (urng == NULL)
    urng = prng_new("mt19937(8888)");

  /* allocate memory for list */
  if (list_of_distr != NULL) {
    *return_list = list_of_distr;
    return n_distr;
  }

  /* else we have to work */
    
  list_of_distr = malloc( N_DISTRIBUTIONS * sizeof(struct list_distr *));

  /* reset counter for distributions */
  n_distr = 0;
  list = list_of_distr;
  
  /**************************************************************************/
  /**                                                                      **/
  /**                Continuous Univariate Distributions                   **/
  /**                                                                      **/
  /**************************************************************************/

#ifdef D_BETA
  /** Beta distribution **/
  fpar[0] = 1.;
  fpar[1] = 2.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 5.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 100.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 100.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1000.;
  fpar[1] = 1000.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 10.;
  fpar[2] = -3.;
  fpar[3] = 15.;
  list->distr = unur_distr_beta(fpar,4);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

#ifdef D_CAUCHY
  /** Cauchy distribution **/
  list->distr = unur_distr_cauchy(NULL,0);
  list->type  = T_Tconcave;
  list->c_max = -0.5;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 20.;
  list->distr = unur_distr_cauchy(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = -0.5;
  ++n_distr; ++list;
#endif

#ifdef D_GAMMA
  /** Gamma distributions **/
  fpar[0] = 1.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 2.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 3.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 10.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1000.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 1.e+10;
  list->distr = unur_distr_gamma(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 1.e-10;
  list->distr = unur_distr_gamma(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 10.;
  fpar[2] = 1.e+10;
  list->distr = unur_distr_gamma(fpar,3);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

#ifdef D_NORMAL
  /** Normal distributions **/
  list->distr = unur_distr_normal(NULL,0);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 10.;
  fpar[1] = 1.e-10;
  list->distr = unur_distr_normal(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 0.;
  fpar[1] = 1.e+10;
  list->distr = unur_distr_normal(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

#ifdef D_UNIFORM
  /** Uniform distribution **/
  list->distr = unur_distr_uniform(NULL,0);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 20.;
  list->distr = unur_distr_uniform(fpar,2);
  list->type  = T_Tconcave;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

  /**************************************************************************/
  /**                                                                      **/
  /**                 Discrete Univariate Distributions                    **/
  /**                                                                      **/
  /**************************************************************************/

#ifdef D_PV_RANDOM
  /** random vector with random entries **/

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_random_vector(10);
  unur_distr_discr_set_prob(distr,prob,10);
  free(prob);
  unur_distr_set_name(distr,"pv(random)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_random_vector(100);
  unur_distr_discr_set_prob(distr,prob,100);
  free(prob);
  unur_distr_set_name(distr,"pv(random)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_random_vector(1000);
  unur_distr_discr_set_prob(distr,prob,1000);
  free(prob);
  unur_distr_set_name(distr,"pv(random)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_random_vector(10000);
  unur_distr_discr_set_prob(distr,prob,10000);
  free(prob);
  unur_distr_set_name(distr,"pv(random)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;
#endif

#ifdef D_PV_GEOMETRIC
  /** random vector with geometrically distributed entries **/

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_geometric_vector(1.,1000);
  unur_distr_discr_set_prob(distr,prob,1000);
  free(prob);
  unur_distr_set_name(distr,"pv(geometric)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_geometric_vector(0.99,1000);
  unur_distr_discr_set_prob(distr,prob,1000);
  free(prob);
  unur_distr_set_name(distr,"pv(geometric)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_geometric_vector(0.9,1000);
  unur_distr_discr_set_prob(distr,prob,1000);
  free(prob);
  unur_distr_set_name(distr,"pv(geometric)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_geometric_vector(0.5,1000);
  unur_distr_discr_set_prob(distr,prob,1000);
  free(prob);
  unur_distr_set_name(distr,"pv(geometric)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

  distr = unur_distr_new(UNUR_DISTR_DISCR);
  prob = make_geometric_vector(0.1,1000);
  unur_distr_discr_set_prob(distr,prob,1000);
  free(prob);
  unur_distr_set_name(distr,"pv(geometric)");
  list->distr = distr;
  list->type  = T_fpv;
  ++n_distr; ++list;

#endif

  /**************************************************************************/
  /**                                                                      **/
  /**                                                                      **/
  /**                                                                      **/
  /**************************************************************************/

  /* check N_DISTRIBUTIONS (compile time constant !!) */
  if (n_distr > N_DISTRIBUTIONS) {
    /* this is very bad */
    fprintf(stderr,"N_DISTRIBUTIONS must be raised to %d\n",n_distr);
    abort();
  }

  *return_list = list_of_distr;
  return n_distr;

} /* end of make_list_of_distributions() */

/*---------------------------------------------------------------------------*/

/* make a probability vector (need not sum to one)
   (use random entries)                                                      */
static double *make_random_vector(int len)
{
  double *prob;
  int i;

  /* allocate memory */
  prob = malloc(len*sizeof(double));
  if (!prob) abort();
     
  /* main part of geometric distribution */
  for( i=0; i<len; i++ ) 
    prob[i] = prng_get_next(urng);

  return prob;

} /* end of make_random_vector() */

/*---------------------------------------------------------------------------*/

/* make a probability vector (need not sum to one)
   (use geometric distribution)                                              */
static double *make_geometric_vector(double q, int len)
{
  double *prob;
  int i;

  /* allocate memory */
  prob = malloc(len * sizeof(double));
  if (!prob) abort();

  /* main part of geometric distribution */
  prob[0] = 1.;
  for( i=1; i<len; i++ ) 
    prob[i] = prob[i-1] * q;

  return prob;

} /* end of make_geometric_vector() */

/*---------------------------------------------------------------------------*/

