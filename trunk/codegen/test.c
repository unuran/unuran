#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unuran.h>

int unurgen( struct unur_gen *gen, FILE *out, const char *distr_name );


int main (void)
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  distr = unur_distr_normal(NULL,0);
  par = unur_tdr_new(distr);
  unur_tdr_set_guidefactor(par,1.);
/*    unur_tdr_set_c(par,0.); */
  gen = unur_init( par );

/*    _unur_tdr_ps_codegen(gen,stdout,NULL); */
  unurgen(gen,stdout,NULL);

  unur_distr_free(distr);
  unur_free(gen);

  exit (0);
}
