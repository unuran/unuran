

#include <unistd.h>
#include <string.h>

#include <unuran.h>
#include <unuran_acg.h>

/* --------------------------------------------------------------------------*/
/* Error codes                                                               */

#define UNUR_ACG_EXIT_SUCCESS       0      /* terminated successfully        */
#define UNUR_ACG_EXIT_USAGE        -1      /* wrong usage                    */
#define UNUR_ACG_EXIT_LANGUAGE      1      /* unknown programming language   */
#define UNUR_ACG_EXIT_DISTRIBUTION  2      /* unknown distribution           */
#define UNUR_ACG_EXIT_DOMAIN        3      /* invalid domain                 */
#define UNUR_ACG_EXIT_PAR          10      /* cannot make parameter object   */
#define UNUR_ACG_EXIT_GEN          20      /* cannot make generator object   */
#define UNUR_ACG_EXIT_CODE         30      /* cannot generate code           */

/* --------------------------------------------------------------------------*/

static void usage (void);
/* --------------------------------------------------------------------------*/
/* Print usage message.                                                      */
/* --------------------------------------------------------------------------*/

int (*get_language(const char *language)) ();
/*---------------------------------------------------------------------------*/
/* Get pointer for code generator for given language.                        */
/* --------------------------------------------------------------------------*/

UNUR_DISTR *(*get_distribution(const char *distribution)) ();
/* --------------------------------------------------------------------------*/
/* Get pointer to make given distribution.                                   */
/* --------------------------------------------------------------------------*/

int get_fparams(char *params, double fparam[]);
/* --------------------------------------------------------------------------*/
/* Get parameters for PDF of distribution.                                   */
/* --------------------------------------------------------------------------*/

int get_domain(char *domain, double fdomain[]);
/* --------------------------------------------------------------------------*/
/* Get domain for distribution.                                              */
/* --------------------------------------------------------------------------*/

/* --------------------------------------------------------------------------*/

static const char progname[] = "acg";
/* --------------------------------------------------------------------------*/
/* program name.                                                             */
/* --------------------------------------------------------------------------*/


/* --------------------------------------------------------------------------*/


/*****************************************************************************/
/**                                                                         **/
/** Main                                                                    **/
/**                                                                         **/
/*****************************************************************************/


int main (int argc, char *argv[]){


  /* pointer to function returning poiter to UNUR_DISTR */
  UNUR_DISTR * (*distrfunc)() = NULL;

  /* pointer to acg, use C as default language */
  int (*langfunc)() = unur_acg_C;

  /* distribution, parameter and generator object */
  UNUR_DISTR *distr;  /* distribution object */
  UNUR_PAR *par;      /* parameter object    */
  UNUR_GEN *gen;      /* generator object    */
  
  /* parameters for distribution */
  double fpar[UNUR_DISTR_MAXPARAMS];
  double domain[2];   /* domain for distribution */

  int domainset = 0; /* checks if domain is set        */

  char c;

  int n_params = 0; /* number of parameters */


  /* ------------------------------------------------------------------------*/
  /* read options                                                            */

  while ((c = getopt(argc, argv, "d:p:l:D:")) != -1) {
    switch (c) {
    case 'd':     /* distribution */
      distrfunc = get_distribution(optarg);
      break;
    case 'p':     /* parameters for PDF */
      n_params = get_fparams(optarg,fpar);
      break;
    case 'D':     /* domain */
      get_domain(optarg,domain);
      domainset = 1;
      break;
    case 'l':     /* progamming language */
      langfunc = get_language(optarg);
      break;
    case '?':    /* Help Message  */
    default:
      usage();
      exit (UNUR_ACG_EXIT_USAGE);
    }
  }

  /* ------------------------------------------------------------------------*/
  /* check whether distribution is set */

  if (distrfunc == NULL){
    fprintf(stderr, "Distribution not set.\n");
    usage();
    exit (UNUR_ACG_EXIT_USAGE);
  }

  /* ------------------------------------------------------------------------*/
  /* create distribution object */

  distr = distrfunc(fpar, n_params); 
  if (distr == NULL) {
    fprintf(stderr, "Cannot create distribution object.\n");
    exit (UNUR_ACG_EXIT_DISTRIBUTION);
  }

  if (domainset == 1) {
    if (!unur_distr_cont_set_domain(distr, domain[0], domain[1])) {
      fprintf(stderr, "Cannot set domain for distribution object.\n");
      exit (UNUR_ACG_EXIT_DOMAIN);
    }
  }

  /* ------------------------------------------------------------------------*/
  /* create generator object                                                 */

  par = unur_tdr_new(distr);       /* parameter object    */
  if (par == NULL) {
    fprintf(stderr, "Cannot create parameter object.\n");
    exit (UNUR_ACG_EXIT_PAR);
  }
  
  /* test version only */
  unur_tdr_set_cpoints(par,4,NULL);

  gen = unur_init( par );          /* generator object    */
  if (gen == NULL) {
    fprintf(stderr, "Cannot create generator object.\n");
    exit (UNUR_ACG_EXIT_GEN);
  }

  /* ------------------------------------------------------------------------*/
  /* generate code                                                           */

  if (!langfunc( gen, stdout, NULL )) {
    fprintf(stderr, "Cannot generate program code.\n");
    exit (UNUR_ACG_EXIT_CODE);
  }

  /* ------------------------------------------------------------------------*/
  /* clear memory and exit                                                   */

  unur_distr_free(distr);
  unur_free(gen);

  exit (UNUR_ACG_EXIT_SUCCESS);

} /* end of main() */


/*****************************************************************************/
/**                                                                         **/
/** Subroutines                                                             **/
/**                                                                         **/
/*****************************************************************************/

void
usage (void)
     /*----------------------------------------------------------------------*/
     /* print usage message.                                                 */
     /*                                                                      */
     /* parameters: none                                                     */
     /*----------------------------------------------------------------------*/
{
  fprintf(stderr,"\n");
  fprintf(stderr,"Usage: %s",progname);
  fprintf(stderr," -d Distribution");
  fprintf(stderr," [-l Language]");
  fprintf(stderr," [-p PDF parameters]");
  fprintf(stderr," [-D Domain]");
  fprintf(stderr,"\n\n");
  fprintf(stderr,"Default for language: C\n");
  fprintf(stderr,"PDF parameters are required for some distributions.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Example:\n");
  fprintf(stderr," To get JAVA code for a generator of the truncatged gamma distribution\n");
  fprintf(stderr," with parameter 3.5 and 4 and domain (0,7) use:\n");
  fprintf(stderr,"\t%s -d beta -p \"3.5 4\" -l Java -D \"0 7\"\n\n",progname);
  fprintf(stderr,"Note parameters and the bounds of the domain are provided\n");
  fprintf(stderr,"within a string seperated by blanks.\n\n");

  exit (UNUR_ACG_EXIT_USAGE);

} /* end of usage() */

/*---------------------------------------------------------------------------*/

int
(*get_language(const char *language)) ()
     /*----------------------------------------------------------------------*/
     /* get pointer for code generator for given language                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*    language ... pointer to string with language description          */
     /*                                                                      */
     /* error:                                                               */
     /*    abort program                                                     */
     /*----------------------------------------------------------------------*/
{
  if ( strcasecmp(language, "JAVA") == 0 ) {
    return unur_acg_JAVA;
  }
  else if ( strcasecmp(language, "C") == 0) {
    return unur_acg_C;
  }
  else if ( strcasecmp(language, "FORTRAN") == 0) {
    return unur_acg_FORTRAN;
  }
  else {
    fprintf(stderr, "Unknown programming language: %s\n",language);
    exit (UNUR_ACG_EXIT_LANGUAGE);
  }

} /* end of get_language() */

/* --------------------------------------------------------------------------*/

UNUR_DISTR *
(*get_distribution(const char *distribution)) ()
     /*----------------------------------------------------------------------*/
     /* get pointer to make given distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*    distribution ... pointer to string with distribution description  */
     /*                                                                      */
     /* error:                                                               */
     /*    abort program                                                     */
     /*----------------------------------------------------------------------*/
{
  if ( strcasecmp(distribution, "normal") == 0 ) {
    return unur_distr_normal;
  }
  else if ( strcasecmp(distribution, "beta") == 0 ) {
    return unur_distr_beta;
  }
  else {
    fprintf(stderr, "Unknown distribution: %s\n",distribution);
    exit (UNUR_ACG_EXIT_DISTRIBUTION);
  }

} /* end of get_distribution() */

/* --------------------------------------------------------------------------*/

int 
get_fparams(char *params, double fparam[])
     /*----------------------------------------------------------------------*/
     /* get parameters for PDF of distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*    params  ... pointer to string with parameters                     */
     /*    fparams ... pointer to array with resulting parameters            */
     /*                                                                      */
     /* return:                                                              */
     /*    number of parameters in list                                      */
     /*                                                                      */
     /* error:                                                               */
     /*    abort program                                                     */
     /*----------------------------------------------------------------------*/
{
  char *toopts;    /* pointer to charakter array */
  char *chktoopts; /* dito, checks correctness   */
  int n_param;

  toopts = params;
  for(n_param=0; n_param<UNUR_DISTR_MAXPARAMS ;n_param++) {
    chktoopts = toopts;
    fparam[n_param] = strtod(toopts, &toopts);

    if (chktoopts == toopts)
      /* no success reading a double */
      break;

  }

  /* o.k. */
  return n_param;

} /* end of get_fparams() */

/* --------------------------------------------------------------------------*/

int 
get_domain(char *domain, double fdomain[])
     /*----------------------------------------------------------------------*/
     /* get domain for PDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*    domain  ... pointer to string with domain                         */
     /*    fdomain ... pointer to array with resulting domain                */
     /*                                                                      */
     /* return:                                                              */
     /*    number of parameters in list                                      */
     /*                                                                      */
     /* error:                                                               */
     /*    abort program                                                     */
     /*----------------------------------------------------------------------*/
{
  char *toopts;    /* pointer to charakter array */
  char *chktoopts; /* dito, checks correctness   */

  chktoopts = toopts = domain;
  fdomain[0] = strtod(toopts, &toopts);
  fdomain[1] = strtod(toopts, &toopts);

  if (chktoopts == toopts || fdomain[0] >= fdomain[1]) {
    fprintf(stderr, "Passed Domain not correct\n");
    exit (UNUR_ACG_EXIT_DOMAIN);
  }

  /* o.k. */
  return 1;

} /* end of get_domain() */

/* --------------------------------------------------------------------------*/

