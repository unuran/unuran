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

#include "t_unuran.h"

/*---------------------------------------------------------------------------*/

struct list_distr *list_of_distr;  /* pointer to list of distributions      */
int n_distr;                        /* number of distributions               */

#define N_DISTRIBUTIONS 100         /* maximum number of distributions       */

/*---------------------------------------------------------------------------*/

/* make list of distributions                                                */
void make_list_of_distributions( void )
{
  double fpar[10];
  struct list_distr *list;

  /* allocate memory for list */
  list_of_distr = malloc( N_DISTRIBUTIONS * sizeof(struct list_distr *));

  /* reset counter for distributions */
  n_distr = 0;
  list = list_of_distr;

#ifdef D_BETA
  /** Beta distribution **/
  fpar[0] = 1.;
  fpar[1] = 2.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 5.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 100.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 100.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1000.;
  fpar[1] = 1000.;
  list->distr = unur_distr_beta(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 10.;
  fpar[2] = -3.;
  fpar[3] = 15.;
  list->distr = unur_distr_beta(fpar,4);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

#ifdef D_CAUCHY
  /** Cauchy distribution **/
  list->distr = unur_distr_cauchy(NULL,0);
  list->type  = T_TYPE_TDR;
  list->c_max = -0.5;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 20.;
  list->distr = unur_distr_cauchy(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = -0.5;
  ++n_distr; ++list;
#endif

#ifdef D_GAMMA
  /** Gamma distributions **/
  fpar[0] = 1.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 2.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 3.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 10.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1000.;
  list->distr = unur_distr_gamma(fpar,1);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 1.e+10;
  list->distr = unur_distr_gamma(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 1.e-10;
  list->distr = unur_distr_gamma(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 5.;
  fpar[1] = 10.;
  fpar[2] = 1.e+10;
  list->distr = unur_distr_gamma(fpar,3);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

#ifdef D_NORMAL
  /** Normal distributions **/
  list->distr = unur_distr_normal(NULL,0);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 10.;
  fpar[1] = 1.e-10;
  list->distr = unur_distr_normal(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 0.;
  fpar[1] = 1.e+10;
  list->distr = unur_distr_normal(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

#ifdef D_UNIFORM
  /** Uniform distribution **/
  list->distr = unur_distr_uniform(NULL,0);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;

  fpar[0] = 1.;
  fpar[1] = 20.;
  list->distr = unur_distr_uniform(fpar,2);
  list->type  = T_TYPE_TDR;
  list->c_max = 0.;
  ++n_distr; ++list;
#endif

  /* check N_DISTRIBUTIONS (compile time constant !!) */
  if (n_distr > N_DISTRIBUTIONS) {
    /* this is very bad */
    fprintf(stderr,"N_DISTRIBUTIONS must be raised to %d\n",n_distr);
    abort();
  }

} /* end of make_list_of_distributions() */

/*---------------------------------------------------------------------------*/

