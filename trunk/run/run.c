/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Examples                                                                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <unuran.h>
#include <unuran_tests.h>

#define RUN_TESTS       (~0x0u)
/* #define RUN_TESTS       UNUR_TEST_SAMPLE */

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr[4];
  UNUR_PAR *par;
  UNUR_GEN *gen;

  UNUR_DISTR *normal;

  int i;
  double p[5], dir[5]; 
  unur_set_default_debug(~0U);

  normal = unur_distr_multinormal(3,NULL,NULL);
  for(i=0;i<3;i++) p[i]=3.*unur_urng_sample(NULL);  
  distr[0] = unur_distr_condi_new( normal, p, NULL, 0 );
  unur_distr_set_name(distr[0],"condi_standardmultinormal_3");
  unur_distr_cont_get_mode(distr[0]);
  unur_distr_free(normal);

  normal = unur_distr_multinormal(3,NULL,NULL);
  for(i=0;i<3;i++) p[i]=3.*unur_urng_sample(NULL);  
  distr[1] = unur_distr_condi_new( normal, p, NULL, 2 );
  unur_distr_set_name(distr[1],"condi_standardmultinormal_3");
  unur_distr_cont_get_mode(distr[1]);
  unur_distr_free(normal);

  normal = unur_distr_multinormal(4,NULL,NULL);
  for(i=0;i<4;i++) p[i]=3.*unur_urng_sample(NULL);  
  for(i=0;i<4;i++) dir[i]=0.5+unur_urng_sample(NULL);  
  distr[2] = unur_distr_condi_new( normal, p, dir, 0 );
  unur_distr_set_name(distr[2],"condi_standardmultinormal_4");
  unur_distr_cont_get_mode(distr[2]);
  unur_distr_free(normal);

  {
    double mean[3], covar[3*3];
    UNUR_DISTR *covar_distr;
    UNUR_GEN *covar_gen;
    UNUR_GEN *mean_gen;
    for(i=0;i<3;i++) p[i]=3.*unur_urng_sample(NULL);  
    for(i=0;i<3;i++) dir[i]=0.5+unur_urng_sample(NULL);  
    mean_gen = unur_str2gen("normal(0,3)");
    for (i=0; i<3; i++) mean[i] = unur_sample_cont(mean_gen);
    unur_free(mean_gen); 
    covar_distr = unur_distr_correlation(3);
    covar_gen = unur_init(unur_mcorr_new(covar_distr));
    do { unur_sample_matr(covar_gen,covar); 
    normal = unur_distr_multinormal(3,mean,covar); 
    } while (normal==NULL);
    unur_distr_free(covar_distr);
    unur_free(covar_gen);
    distr[3] = unur_distr_condi_new( normal, p, dir, 0 );
    unur_distr_set_name(distr[3],"condi_multinormal_random");
    unur_distr_cont_get_mode(distr[3]);
    unur_distr_free(normal);
  }

  for (i=0; i<4;i++) {
    par = unur_tdr_new(distr[i]);
/*     gen = unur_init(par); */
    unur_test_par_count_pdf( par, 100000, 2, stdout );
/*     unur_free(gen); */
  }

  return 0;
}

/*---------------------------------------------------------------------------*/

