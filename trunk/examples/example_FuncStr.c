/* ------------------------------------------------------------- */
/* File: example_FuncStr.c                                       */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

/* ------------------------------------------------------------- */

/* Example how to define the PDF for a continuous univariate     */
/* distribution using a function string.                         */

/* ------------------------------------------------------------- */

int main()
{
  UNUR_DISTR *distr;    /* distribution object */

  /* Get empty distribution object for a continuous distribution */
  distr = unur_distr_cont_new();

  /* Set PDF using function string */
  unur_distr_cont_set_pdfstr(distr,"1-x*x");
  unur_distr_cont_set_domain(distr,-1.,1.);
 
  /* Read function string */
  printf("functionstring: %s\n",unur_distr_cont_get_pdfstr(distr));

  /* Destroy distribution object */
  unur_distr_free(distr);

  exit (EXIT_SUCCESS);

} /* end of main() */

/* ------------------------------------------------------------- */
