#!/usr/bin/perl
# ----------------------------------------------------------------
# make C programm to run Automatic Code Generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

require "read_PDF.pl";

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata('..');

# For description of data fields in this list see file `readPDF.pl'.

# ----------------------------------------------------------------
# Print C source

print acg_header();
print acg_main();
print acg_usage();
print acg_distribution();
print acg_fparams();
print acg_domain();
print acg_n_cpoints();
print acg_language();
print acg_logfile();

# ----------------------------------------------------------------
exit 0;
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# make header for C file 
# ----------------------------------------------------------------

sub acg_header
{
    my $header = <<EOS;

/* ------------------------------------------------------------------------- */
/* Include header file                                                       */
/* ------------------------------------------------------------------------- */

\#include <unistd.h>
\#include <string.h>

\#include <unuran.h>
\#include <unuran_acg.h>

/* ------------------------------------------------------------------------- */
/* Program name                                                              */
/* ------------------------------------------------------------------------- */

static const char progname[] = "acg";

/* ------------------------------------------------------------------------- */
/* Exit codes                                                                */
/* ------------------------------------------------------------------------- */

\#define ACG_EXIT_SUCCESS       0     /* terminated successfully              */
\#define ACG_EXIT_FAIL_INPUT   -1     /* input error                          */
\#define ACG_EXIT_FAIL_GEN      1     /* cannot create generator object       */
\#define ACG_EXIT_FAIL_CODE     2     /* cannot generator code                */

/* ------------------------------------------------------------------------- */
/* Local prototypes                                                          */
/* ------------------------------------------------------------------------- */

/* Print usage message.                                                      */
void usage (void);

/* Get pointer for code generator for given language.                        */
int (*get_language(const char *language)) ();

/* Get pointer to make given distribution.                                   */
UNUR_DISTR *(*get_distribution(const char *distribution)) ();

/* Get parameters for PDF of distribution.                                   */
int get_fparams(char *params, double fparam[]);

/* Get domain for distribution.                                              */
int get_domain(char *domain, double fdomain[]);

/* Get number of construction points.                                        */
int get_n_cpoints(char *n_cpoints);

/* Get name of log file                                                      */
const char *get_logfile(const char *logfile);

/* ------------------------------------------------------------------------- */
/* Error messages                                                            */
/* ------------------------------------------------------------------------- */

EOS

    if ($DEBUG) {
	$header .= "#define fatal(msg) do { fprintf(stderr,msg); } while(0)\n\n";
    }
    else {
	$header .= "#define fatal(msg) do { } while(0)\n\n";
    }

    return $header;
} # end of acg_header()

# ----------------------------------------------------------------
# make main
# ----------------------------------------------------------------

sub acg_main {

    my $main = <<EOS;

/*****************************************************************************/
/**                                                                         **/
/** Main                                                                    **/
/**                                                                         **/
/*****************************************************************************/


int main (int argc, char *argv[]){

  /* pointer to logfile */
  const char *logfile = NULL;
  FILE *logstream = NULL;

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

  int domainset = 0;  /* checks if domain is set */
  int n_params = 0;   /* number of parameters */
  int n_cpoints = 30; /* number of construction points (default value) */

  char c;

  /* ------------------------------------------------------------------------*/
  /* read options                                                            */

  while ((c = getopt(argc, argv, "d:p:D:n:l:L:")) != -1) {
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
    case 'n':     /* number of construction points */
      n_cpoints = get_n_cpoints(optarg);
      break;
    case 'l':     /* progamming language */
      langfunc = get_language(optarg);
      break;
    case 'L':
      logfile = get_logfile(optarg);
      break;
    case '?':    /* Help Message  */
    default:
      usage();
      exit (ACG_EXIT_FAIL_INPUT);
    }
  }

  /* ------------------------------------------------------------------------*/
  /* check whether distribution is set */

  if (distrfunc == NULL) {
    fatal("Distribution not set.\\n");
    usage();
    exit (ACG_EXIT_FAIL_INPUT);
  }

  /* ------------------------------------------------------------------------*/
  /* logging */

  if (logfile) {
      logstream = fopen(logfile,"w");
      unur_set_stream(logstream);
      unur_set_default_debug(UNUR_DEBUG_ALL); 
  }
  else   
      unur_set_default_debug(UNUR_DEBUG_OFF); 

  /* ------------------------------------------------------------------------*/
  /* create distribution object */

  distr = distrfunc(fpar, n_params); 
  if (distr == NULL) {
    fatal("Cannot create distribution object.\\n");
    exit (ACG_EXIT_FAIL_INPUT);
  }

  if (domainset == 1) {
    if (!unur_distr_cont_set_domain(distr, domain[0], domain[1])) {
      fatal("Cannot set domain for distribution object.\\n");
      exit (ACG_EXIT_FAIL_INPUT);
    }
  }

  /* ------------------------------------------------------------------------*/
  /* create generator object                                                 */

  par = unur_tdr_new(distr);       /* parameter object    */
  if (par == NULL) {
    fatal("Cannot create parameter object.\\n");
    exit (ACG_EXIT_FAIL_INPUT);
  }
  
  /* test version only */
  unur_tdr_set_cpoints(par,n_cpoints,NULL);

  gen = unur_init( par );          /* generator object    */
  if (gen == NULL) {
    fatal("Cannot create generator object.\\n");
    exit (ACG_EXIT_FAIL_GEN);
  }

  /* ------------------------------------------------------------------------*/
  /* generate code                                                           */

  if (!langfunc( gen, stdout, NULL )) {
    fatal("Cannot generate program code.\\n");
    exit (ACG_EXIT_FAIL_CODE);
  }

  /* ------------------------------------------------------------------------*/
  /* clear memory and exit                                                   */

  unur_distr_free(distr);
  unur_free(gen);

  exit (ACG_EXIT_SUCCESS);

} /* end of main() */


/*****************************************************************************/
/**                                                                         **/
/** Subroutines                                                             **/
/**                                                                         **/
/*****************************************************************************/


EOS

    return $main;
} # end of acg_main()

# ----------------------------------------------------------------
# usage
# ----------------------------------------------------------------

sub acg_usage {

    my $usage = <<EOS;

void
usage (void)
     /*----------------------------------------------------------------------*/
     /* print usage message.                                                 */
     /*                                                                      */
     /* parameters: none                                                     */
     /*----------------------------------------------------------------------*/
{
  fprintf(stderr,"\\n");
  fprintf(stderr,"Usage: %s",progname);
  fprintf(stderr," -d Distribution");
  fprintf(stderr," [-p PDF parameters]");
  fprintf(stderr," [-D Domain]");
  fprintf(stderr," [-n number of construction points]");
  fprintf(stderr," [-l Language]");
  fprintf(stderr," [-L log file]");
  fprintf(stderr,"\\n\\n");
  fprintf(stderr,"Default for language: C\\n");
  fprintf(stderr,"PDF parameters are required for some distributions.\\n");
  fprintf(stderr,"\\n");
  fprintf(stderr,"Example:\\n");
  fprintf(stderr," To get JAVA code for a generator of the truncatged gamma distribution\\n");
  fprintf(stderr," with parameter 3.5 and 4 and domain (0,0.4) use:\\n");
  fprintf(stderr,"\\t%s -d beta -p \\\"3.5 4\\\" -l Java -D \\\"0 0.4\\\"\\n\\n",progname);
  fprintf(stderr,"Note parameters and the bounds of the domain are provided\\n");
  fprintf(stderr,"within a string seperated by blanks.\\n\\n");

  exit (ACG_EXIT_FAIL_INPUT);

} /* end of usage() */

/*---------------------------------------------------------------------------*/

EOS

    return $usage;
} # end of acg_usage()

# ----------------------------------------------------------------
# distribution
# ----------------------------------------------------------------

sub acg_distribution {

    my $distribution = <<EOS;

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
EOS

    # list of distributions
    foreach my $d (sort keys %{$DISTR}) { 
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	$distribution .= "  if ( strcasecmp(distribution, \"$d\") == 0 ) {\n";
	$distribution .= "    return unur_distr_$d;\n  }\n";
    }


    $distribution .= <<EOS;
  else {
    fatal("Unknown distribution.\\n");
    exit (ACG_EXIT_FAIL_INPUT);
  }

} /* end of get_distribution() */

/* --------------------------------------------------------------------------*/

EOS

    return $distribution;
} # acg_distribution()

# ----------------------------------------------------------------
# PDF parameters
# ----------------------------------------------------------------

sub acg_fparams {

    my $fparams = <<EOS;

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

EOS

    return $fparams;
} # end of acg_fparams()

# ----------------------------------------------------------------
# domain
# ----------------------------------------------------------------

sub acg_domain {

    my $domain = <<EOS;

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
    fatal("Passed Domain not correct\\n");
    exit (ACG_EXIT_FAIL_INPUT);
  }

  /* o.k. */
  return 1;

} /* end of get_domain() */

/* --------------------------------------------------------------------------*/

EOS

    return $domain;
} # end of acg_domain()

# ----------------------------------------------------------------
# number of construnction points 
# ----------------------------------------------------------------

sub acg_n_cpoints {

    my $n_cpoints = <<EOS;

int
get_n_cpoints(char *n_cpoints)
     /*----------------------------------------------------------------------*/
     /* Get number of construction points.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*    n_cpoints ... pointer to string with number of cpoints            */
     /*                                                                      */
     /* return:                                                              */
     /*    number of construction points                                     */
     /*                                                                      */
     /* error:                                                               */
     /*    abort program                                                     */
     /*----------------------------------------------------------------------*/
{
  char *toopts;    /* pointer to charakter array */
  char *chktoopts; /* dito, checks correctness   */
  int n_cp;

  chktoopts = toopts = n_cpoints;
  n_cp = strtol(toopts, &toopts,10);

  if (chktoopts == toopts || n_cp < 3) {
    fatal("Passed n_cpoints not correct, use default instead.\\n");
    n_cp = 30;
  }

  /* o.k. */
  return n_cp;

} /* end of get_n_cpoints() */

/* --------------------------------------------------------------------------*/

EOS

    return $n_cpoints;
} # end of acg_n_cpoints()

# ----------------------------------------------------------------
# parse language
# ----------------------------------------------------------------

sub acg_language {

    my $language = <<EOS;

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
    fatal("Unknown programming language.\\n");
    exit (ACG_EXIT_FAIL_INPUT);
  }

} /* end of get_language() */

/* --------------------------------------------------------------------------*/

EOS

    return $language;
} # end of acg_language()

# ----------------------------------------------------------------
# parse log file name
# ----------------------------------------------------------------

sub acg_logfile {

    my $logfile = <<EOS;

const char *
get_logfile(const char *logfile)
     /*----------------------------------------------------------------------*/
     /* get pointer to log file name                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*    logfile ... pointer to string with log file name                  */
     /*                                                                      */
     /* error:                                                               */
     /*    abort program                                                     */
     /*----------------------------------------------------------------------*/
{
    return logfile;
} /* end of get_language() */

/* --------------------------------------------------------------------------*/

EOS

    return $logfile;
} # end of acg_logfile()

# ----------------------------------------------------------------
# end
# ----------------------------------------------------------------









