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

  fpar[0] = 2.;
  list->distr = unur_distr_student(fpar,1);
  list->type  = T_TYPE_TDR;
  list->c_max = -0.5;
  ++n_distr; ++list;


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

  /* check N_DISTRIBUTIONS (compile time constant !!) */
  if (n_distr > N_DISTRIBUTIONS) {
    /* this is very bad */
    fprintf(stderr,"N_DISTRIBUTIONS must be raised to %d\n",n_distr);
    abort();
  }

} /* end of make_list_of_distributions() */

/*---------------------------------------------------------------------------*/

