#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unuran.h>


extern void         subst_whitespace(char *, char);
extern char       * elim_whitespace(char *);
extern int          parselist(char *, double *);
extern UNUR_DISTR * make_distr_obj(char *);
extern UNUR_PAR   * make_par_obj(UNUR_DISTR *, char *);
extern UNUR_GEN   * make_gen_obj (char *);

/* main -- test program */
int main(){

  int i;

  UNUR_GEN *gen;
  char str[] = "distr=normal(-1, .5); domain=(-1,-.5) : method=tdr; c=-0.5; variant_ia";
  gen = make_gen_obj(str);

  for ( i=0; i<15; i++){
    printf("rand num: %f\n", unur_sample_cont(gen));
  }

  return (0);
}
