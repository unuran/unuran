
/* ------------------------------------------------------------------------- */
/* Include header file                                                       */
/* ------------------------------------------------------------------------- */

#include <unistd.h>
#include <string.h>

#include <unuran.h>

/* ------------------------------------------------------------------------- */
/* Program name                                                              */
/* ------------------------------------------------------------------------- */

static const char progname[] = "rungen";

/* ------------------------------------------------------------------------- */
/* Local prototypes                                                          */
/* ------------------------------------------------------------------------- */

/* Print usage message.                                                      */
void usage (void);

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

int main (int argc, char *argv[])
{
  /* string with description of generator for String API */
  const char *gen_str = NULL;
  
  /* number points generated */
  int n = 0;
  
  /* pointer to logfile */
  const char *logfile = NULL;
  FILE *logstream = NULL;

  /* generator object */
  UNUR_GEN *gen;       /* generator object    */
  
  char c;

  /* ------------------------------------------------------------------------*/
  /* read options                                                            */

  while ((c = getopt(argc, argv, "L:n:")) != -1) {
    switch (c) {
    case 'L':     /* name of log file */
      logfile = optarg;
      break;
    case 'n':
      n = atoi(optarg);
      break;
    case '?':     /* Help Message  */
    default:
      usage();
      exit (EXIT_FAILURE);
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
    exit (EXIT_FAILURE);
  }

  /* ------------------------------------------------------------------------*/
  /* logging                                                                 */

  if (logfile) {
    /* open log file */
    if (strcmp("stdout",logfile))
      logstream = fopen(logfile,"w");
    else /* use stdout instead of a log file */
      logstream = stdout;

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
    exit (EXIT_FAILURE);
  }

  /* ------------------------------------------------------------------------*/
  /* run generator                                                           */
  while (_unur_tdr_is_ARS_running(gen))
    unur_sample_cont(gen);

  /* ------------------------------------------------------------------------*/
  /* clear memory and exit                                                   */

  unur_free(gen);

  exit (EXIT_SUCCESS);

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
  fprintf(stderr," [-L log file]");
  fprintf(stderr," \"generator string\"");
  fprintf(stderr,"\n\n");
  fprintf(stderr," -L ... log file for debugging information.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"For syntax of generator string see UNURAN Manual, String API.\n");
} /* end of usage() */

/*---------------------------------------------------------------------------*/
