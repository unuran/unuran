/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/
    
/*---------------------------------------------------------------------------*/

/*  #define DEBUG 1 */

/*---------------------------------------------------------------------------*/

#include "test_with_Mathematica.h"

/*---------------------------------------------------------------------------*/

int main()
{ 
  double fpm[10]={1.,2.,3.};

  UNUR_DISTR *distr;    /* distribution */

  distr = unur_distr_beta(fpm,2);
  test_cont_cdf_pdf( distr, "t_distr_beta.data", 1. );
  unur_distr_free(distr);

  distr = unur_distr_gamma(fpm,2);
  test_cont_cdf_pdf( distr, "t_distr_gamma.data", 1. );
  unur_distr_free(distr);

  distr = unur_distr_normal(fpm,2);
  test_cont_cdf_pdf( distr, "t_distr_normal.data", 1. );
  unur_distr_free(distr);

  exit (0);
}

/*---------------------------------------------------------------------------*/





