#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <unuran.h>
#include <unuran_acg.h>
#include <unuran_tests.h>


extern int strcpy();
extern int strcasecmp();

/* ********************************************************* */
/*                                                           */
/* ********************************************************* */
int main (int argc, char *argv[]){


  /* pointer to function returning poiter to UNUR_DISTR */
  UNUR_DISTR * (*distrfunc)() = NULL;

  /* pointer to function returning nothing */
  int (*langfunc)() = NULL;

  /* distribution, parameter and generator object */
  UNUR_DISTR *distr;  /* distribution object */
  UNUR_PAR *par;      /* parameter object    */
  UNUR_GEN *gen;      /* generator object    */
  
  /* parameters for distribution */
  double fpar[10];  /* max 10 parameters for distribution */
  double domain[2]; /* domain for distribution */

  int distrset = 0;  /* checks if distribution is set  */
  int paramsset = 0; /* checks if parameters are  set  */
  int langset = 0;   /* checks if prog language is set */
  int domainset = 0; /* checks if domain is set        */

  char c;
  char distropt[63];  /* holds the distribution    */
  char paramsopt[63]; /* holds the parameters      */
  char langopt[63];   /* holds the output language */
  char domainopt[63]; /* holds the domain          */

  char *toopts;    /* pointer to charakter array */
  char *chktoopts; /* dito, checks correctness   */

  int no_of_params; /* number of parameters */
  int i;            /* just a counter       */


  /* ----------------------------------------------------------------------*/
  /* read options                                                          */
  int exit_readopts = 0;
  while ((c = getopt(argc, argv, "d:p:l:D:")) != -1) {
    switch (c) {
    case 'd':
      strcpy(distropt, optarg);
      distrset = 1;
       break;
    case 'p':
      strcpy(paramsopt, optarg);
      paramsset = 1;
      break;
    case 'l':
      strcpy(langopt, optarg);
      langset = 1;
      break;
    case 'D':
      strcpy(domainopt, optarg);
      domainset = 1;
      break;
    case '?':    /* Help Message  */
    default:
      printf("\nUsage of thisfunc:\n   thisfunc -l language");
      printf(" -d Distribution [-p parameters] [-D Domain]\n");
      printf("\nExample:\nTo get JAVA code of the beta distribution ");
      printf("with the parameters 5, 3.5 within the domain (-1, 5) use this");
      printf("function call:\n");
      printf("\n   thisfunc -d beta -p \"5 3.5\" -l Java -D \"-1.0 5\"\n\n");
      printf("Note the parameters and the bounds of the domain are provided");
      printf(" within a string seperated by blanks.\n\n");

      exit_readopts = 1;
      return 1;
    }
    if (exit_readopts == 1) break;
  }

  /* check whether distribution and output language are set */
  if (distrset == 0 || langset == 0){
    fprintf(stderr, "Distribution or language not set\n");
    return 1;
  }

  /* end of read options                                                   */
  /* --------------------------------------------------------------------- */



  /* ------------------------------------------------------------- */
  /* determine distribution */
  if ( strcasecmp(distropt, "normal") == 0 ){
    distrfunc = unur_distr_normal;
  }
  else if(strcasecmp(distropt, "beta") == 0){
    distrfunc = unur_distr_beta;
  }
  else{
    fprintf(stderr, "Unknown distribution\n");
  }

  /* ------------------------------------------------------------- */
  /* determine programming language */
  if ( strcasecmp(langopt, "java") == 0 ){
    langfunc = unur_acg_JAVA;
  }
  else if ( strcasecmp(langopt, "c") == 0){
    langfunc = unur_acg_C;
  }
  else if ( strcasecmp(langopt, "fortran") == 0){
    langfunc = unur_acg_FORTRAN;
  }
  else{
    fprintf(stderr, "Unknown programming language\n");
  }

  /* ------------------------------------------------------------- */
  /* parse and set parameters */
  toopts = paramsopt;
  for(i=0; i<10 ;i++){
    chktoopts = toopts;

    fpar[i] = strtod(toopts, &toopts);

    /* no success reading a double */
    if (chktoopts == toopts)
      break;
  }
  no_of_params = i;  /* because of break not incremented */
  /* now fpar contains the parameters and
     no_of_parameters the number of parameters           */

  /* ------------------------------------------------------------- */
  /* parse and set domain */
  toopts = domainopt;
  if (domainset == 1){
    chktoopts = toopts;
    domain[0] = strtod(toopts, &toopts);
    domain[1] = strtod(toopts, &toopts);
    if (chktoopts == toopts){
      fprintf(stderr, "Passed Domain not correct\n");
    }

  }

  /* ------------------------------------------------------------- */
  /* Create Objects and code */

  distr = distrfunc(fpar, no_of_params); /* distribution object */
  par = unur_tdr_new(distr);       /* parameter object    */
  unur_tdr_set_cpoints(par,4,NULL);
  gen = unur_init( par );          /* generator object    */

  /* set domain if specified */
  if (domainset == 1)
    unur_distr_cont_set_domain(distr, domain[0], domain[1]);

  /* language dependent output */
  langfunc( gen, stdout, NULL );

  unur_distr_free(distr);
  unur_free(gen);

  exit (0);

}









