
/* ------------------------------------------------------------------------- */
/* Include header file                                                       */
/* ------------------------------------------------------------------------- */

#include <unistd.h>
#include <string.h>

#include <unuran.h>
#include <unuran_acg.h>

/* ------------------------------------------------------------------------- */
/* Program name                                                              */
/* ------------------------------------------------------------------------- */

static const char progname[] = "acg";

/* ------------------------------------------------------------------------- */
/* Exit codes                                                                */
/* ------------------------------------------------------------------------- */

#define ACG_EXIT_SUCCESS       0     /* terminated successfully              */
#define ACG_EXIT_FAIL_INPUT   -1     /* input error                          */
#define ACG_EXIT_FAIL_GEN      1     /* cannot create generator object       */
#define ACG_EXIT_FAIL_CODE     2     /* cannot generator code                */

/* ------------------------------------------------------------------------- */
/* Local prototypes                                                          */
/* ------------------------------------------------------------------------- */

/* Print usage message.                                                      */
void usage (void);

/* Get pointer for code generator for given language.                        */
int (*get_language(const char *language)) ();

/* ------------------------------------------------------------------------- */
/* Error messages                                                            */
/* ------------------------------------------------------------------------- */

static int debug = 0;
#define fatal(msg) do { if (debug) fprintf(stderr,msg); } while (0)

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/**                                                                         **/
/** Main                                                                    **/
/**                                                                         **/
/*****************************************************************************/

int main (int argc, char *argv[]){

  /* whether to include main into source */
  int with_main = 0;

  /* string with description of generator for String API */
  const char *gen_str = NULL;

  /* name for generator routine */
  const char *gen_name = NULL;

  /* pointer to logfile */
  const char *logfile = NULL;
  FILE *logstream = NULL;

  /* pointer to acg, use C as default language */
  int (*langfunc)() = unur_acg_C;

  /* generator object */
  UNUR_GEN *gen;       /* generator object    */
  
  char c;

  /* ------------------------------------------------------------------------*/
  /* read options                                                            */

  while ((c = getopt(argc, argv, "N:l:L:DM")) != -1) {
    switch (c) {
    case 'N':     /* name for generator routine */
      gen_name = optarg;
      break;
    case 'l':     /* progamming language */
      langfunc = get_language(optarg);
      break;
    case 'L':     /* name of log file */
      logfile = optarg;
      break;
    case 'M':     /* include main */
      with_main = 1;
      break;
    case 'D':     /* print error message on stderr (Debug mode) */
      debug = 1;
      break;
    case '?':     /* Help Message  */
    default:
      usage();
      exit (ACG_EXIT_FAIL_INPUT);
    }
  }

  /* ------------------------------------------------------------------------*/
  /* read generator string                                                   */

  if (optind < argc) {
    gen_str = argv[optind];
  }
  else {
    fatal("No generator string given.\n");
    usage();
    exit (ACG_EXIT_FAIL_INPUT);
  }

  /* ------------------------------------------------------------------------*/
  /* logging                                                                 */

  if (logfile) {
      logstream = fopen(logfile,"w");
      unur_set_stream(logstream);
      unur_set_default_debug(UNUR_DEBUG_ALL); 
  }
  else   
      unur_set_default_debug(UNUR_DEBUG_OFF); 

  /* ------------------------------------------------------------------------*/
  /* create generator object                                                 */
  gen = unur_str2gen(gen_str);         

  if (gen == NULL) {
    fatal("Cannot create generator object.\n");
    exit (ACG_EXIT_FAIL_GEN);
  }

  /* ------------------------------------------------------------------------*/
  /* generate code                                                           */

  if (!langfunc( gen, stdout, gen_name, with_main )) {
    fatal("Cannot generate program code.\n");
    exit (ACG_EXIT_FAIL_CODE);
  }

  /* ------------------------------------------------------------------------*/
  /* clear memory and exit                                                   */

  unur_free(gen);

  exit (ACG_EXIT_SUCCESS);

} /* end of main() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** Subroutines                                                             **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

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
  fprintf(stderr," [-l Language]");
  fprintf(stderr," [-L log file]");
  fprintf(stderr," [-D]");
  fprintf(stderr," [-M]");
  fprintf(stderr," \"generator string\"");
  fprintf(stderr,"\n\n");
  fprintf(stderr," -l ... programming language (C|FORTRAN|JAVA). Default is C.\n");
  fprintf(stderr," -L ... log file for debugging information.\n");
  fprintf(stderr," -D ... debugging mode.\n");
  fprintf(stderr," -M ... inklude main().\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"For syntax of generator string see UNURAN Manual, String API.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Example:\n");
  fprintf(stderr," To get JAVA code for a generator of the truncatged beta distribution\n");
  fprintf(stderr," with parameter 3.5 and 4 and domain (0,0.4) use:\n");
  fprintf(stderr,"\t%s -l Java \"beta(3.5,4); domain=(0,0.4)\"\n\n",progname);
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
  if ( strcasecmp(language, "C") == 0) {
    return unur_acg_C;
  }
  else if ( strcasecmp(language, "FORTRAN") == 0) {
    return unur_acg_FORTRAN;
  }
  else if ( strcasecmp(language, "JAVA") == 0 ) {
    return unur_acg_JAVA;
  }
  else if ( strcasecmp(language, "UNURAN") == 0) {
    return unur_acg_UNURAN;
  }
  else {
    fatal("Unknown programming language.\n");
    exit (ACG_EXIT_FAIL_INPUT);
  }

} /* end of get_language() */

/* --------------------------------------------------------------------------*/
