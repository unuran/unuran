#include <unuran.h>

int main()
{
  UNUR_DISTR *distr;    /* distribution object   */
  char       *functionstring = "1-x*x";

  distr = unur_distr_cont_new();
  unur_distr_cont_set_pdfstr(distr,functionstring); 
 

  printf("functionstring_:\n%s\n",unur_distr_cont_get_pdfstr(distr));
  unur_distr_free(distr);

  exit (EXIT_SUCCESS);
} 
