#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unuran.h>


extern void         subst_whitespace(char *, char);
extern char       * elim_whitespace(char *);
extern int          parselist(char *, double *);
extern UNUR_DISTR * _unur_str2distr(char *);
extern UNUR_PAR   * _unur_str2par(UNUR_DISTR *, char *);
extern UNUR_GEN   * unur_str2gen(char *);

/* main -- test program */
int main(){

  char str[1024];
  int i;

  UNUR_GEN *gen;

  /* test1 */
  strcpy(str, "distr=normal(0, 1):prng = MT19937(123): method=tdr; c=-0.5; variant_ia");
  printf("%s\n", str);
  gen = unur_str2gen(str);
  for ( i=0; i<5; i++){
    printf("rand num: %f\n", unur_sample_cont(gen));
  }
  /* test2 */
  strcpy(str, "distr=normal(0, 1): method=ninv;");
  printf("%s\n", str);
  gen = unur_str2gen(str);
  for ( i=0; i<5; i++){
    printf("rand num: %f\n", unur_sample_cont(gen));
  }
  /* test3 */
  strcpy(str, "distr=normal(0, 1):prng = MT19937(6)");
  printf("%s\n", str);
  gen = unur_str2gen(str);
  for ( i=0; i<5; i++){
    printf("rand num: %f\n", unur_sample_cont(gen));
  }
  /* test4 */
  strcpy(str, "distr=beta(1, 1): method=ninv; c=-0.5; variant_ia:prng=MT19937(454)");
  printf("%s\n", str);
  gen = unur_str2gen(str);
  for ( i=0; i<5; i++){
    printf("rand num: %f\n", unur_sample_cont(gen));
  }

  return (0);
}
