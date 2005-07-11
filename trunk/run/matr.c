#include <unuran.h>
#include <unur_struct.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <utils/vector_source.h>
#include <utils/matrix_source.h>
#include "../tests/testunuran.h"
#include "../tests/testdistributions/testdistributions.h"

/*----------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
#define idx(i,j) ((i)*dim+(j))

  UNUR_DISTR *distr;
  UNUR_PAR   *par;
  UNUR_GEN   *gen;

  int i,j;
  int dim = 2 ;
  double rho = 0.2;
  
  double *mean, *x;
  double *covar; 
  double *covar_inv;
  double *product;
  double *eigenvectors;
  double det;

  if (argc==2) {dim = atoi(argv[1]);}
  if (argc==3) {dim = atoi(argv[1]); rho=atof(argv[2]); }

  x = _unur_vector_new(dim);
  mean = _unur_vector_new(dim);
//  covar = _unur_vector_new(dim*dim);
//  covar_inv = _unur_vector_new(dim*dim);
  product = _unur_vector_new(dim*dim);
  eigenvectors = _unur_vector_new(dim*dim);

 //distr = unur_distr_multinormal_ar1( dim, mean, rho );
 //distr = unur_distr_multinormal_constant_rho( dim, mean, rho );

 distr = unur_distr_multicauchy_ar1( dim, mean, rho );
 
  
/*  
  for (i=0; i<dim; i++)
  for (j=0; j<dim; j++)
    covar[idx(i,j)]=1.;
*/  
 
 covar=unur_distr_cvec_get_covar(distr);
 covar_inv=unur_distr_cvec_get_covar_inv(distr);
 
 _unur_matrix_print_matrix(dim, covar, "COVAR", stdout , "","");
 _unur_matrix_print_matrix(dim, covar_inv, "COVAR INV _unur_distr_cvec_get_covar_inv()", stdout , "",""); 
 _unur_matrix_multiplication(dim, covar, covar_inv, product);
 _unur_matrix_print_matrix(dim, product, "PRODUCT", stdout , "",""); 
 
 memset(covar_inv, 0,dim*dim*sizeof(double));
 _unur_matrix_invert_matrix(dim, covar, covar_inv, &det);
 _unur_matrix_print_matrix(dim, covar_inv, "COVAR INV _unur_matrix_invert_matrix()", stdout , "",""); 
 _unur_matrix_multiplication(dim, covar, covar_inv, product);
 _unur_matrix_print_matrix(dim, product, "PRODUCT", stdout , "",""); 

  printf("determinant=%e\n", det);

//  par=unur_vmt_new(distr);
  par=unur_vnrou_new(distr);
  unur_vnrou_set_r(par, .5);
  gen=unur_init(par);

  if (gen!=NULL) {
    printf("ok\n");
#if 1
//    unur_urng_reset(unur_get_default_urng());
    for (i=1; i<=10; i++) {
      unur_sample_vec(gen, x);
      _unur_matrix_print_vector(dim, x, "",stdout,"","");
    }
#endif
  }

  unur_distr_free(distr);
  unur_free(gen);

  free(eigenvectors);
  free(mean); free(x); 
  free(product);
  //free(covar); free(covar_inv);

  return 0;
}


