/***************************************************************************** 
 *                                                                           * 
 *          UNURAN -- Universal Non-Uniform Random number generator          * 
 *                                                                           * 
 *****************************************************************************/
                                                                               
/*---------------------------------------------------------------------------*/
/*  #define DEBUG 1 */                                                         
/*---------------------------------------------------------------------------*/
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <unuran.h>
#include <config.h>

#ifdef WITH_DMALLOC
#  include <dmalloc.h>
#endif
/*---------------------------------------------------------------------------*/

/* name of data file */
static const char datafile[] = "t_functionparser.data";

/* number of failed tests */
static int n_failed = 0;    

/* file handles */
FILE *DATA;            /* file with data    */
FILE *UNURANLOG;       /* unuran log file   */
FILE *TESTLOG;         /* test log file     */

/*---------------------------------------------------------------------------*/

static int _unur_test_function (void);

/*---------------------------------------------------------------------------*/
                                                                               
int main()                                                                     
{                                                                              
  /* open log file for unuran and set output stream for unuran messages */    
  if ( (UNURANLOG = fopen( "t_functionparser_unuran.log","w" )) == NULL )
    exit (-1);                                           
  unur_set_stream( UNURANLOG );
                                                                               
  /* open log file for testing */                                             
  if ( (TESTLOG = fopen( "t_functionparser_test.log","w" )) == NULL )
    exit (-1);                                             
                                                                               
  /* write header into log file */                                            
  {                                                                           
    time_t started;                                                          
    fprintf(TESTLOG,"\nUNURAN - Universal Non-Uniform RANdom number generator\n\n");
    if (time( &started ) != -1)                                              
      fprintf(TESTLOG,"%s",ctime(&started));                              
    fprintf(TESTLOG,"\n======================================================\n\n");
    /*      fprintf(TESTLOG,"(Search for string \"data\" to find new section.)\n\n"); */
  }                                                                           

  /* open data file */
  if ( (DATA = fopen( datafile,"r" )) == NULL ) {
    printf("ERROR: could not open file %s \n", datafile);
    fprintf(TESTLOG,"ERROR: could not open file %s \n", datafile);
    exit (-1);                                             
  }

  while (_unur_test_function());

  /* close files */
  fclose (UNURANLOG);
  fclose (TESTLOG);
  fclose (DATA);
                               
  /* end */
  exit ((n_failed > 0) ? -1 : 0);
}                                                                              

/*---------------------------------------------------------------------------*/

int
_unur_test_function (void)
{
#define BUFSIZE      1024       /* size of line buffer */

  char buffer[BUFSIZE];         /* buffer for reading line */
  char *ptr_buffer;             /* pointer into buffer */

  char fstr[BUFSIZE];           /* function string */

  double x;                     /* argument x */
  double fx_exp, dfx_exp;       /* expected values for f(x) and f'(x) */
  double fx_obs, dfx_obs;       /* observed values for f(x) and f'(x) */

  UNUR_DISTR *distr;

  /* find next function string */
  while (1) {
    /* read next line */
    fgets(buffer, BUFSIZE, DATA);
    if (feof(DATA)) return 0;
    if (strncmp(buffer,"function=",8)==0) break;
  }

  /* get function string */
  strcpy( fstr, buffer+9 );
  /* remove newline character */
  fstr[strlen(fstr)-1] = '\0';

  /* print info into log file */
  fprintf(TESTLOG,"\"%s\" ...\n", fstr);

  /* make distribution object with given function as PDF */
  distr = unur_distr_cont_new();
  if (unur_distr_cont_set_pdfstr(distr,fstr) == 0) {
    printf("ERROR: syntax error in \"%s\"\n", fstr);
    fprintf(TESTLOG,"ERROR: syntax error in \"%s\"\n", fstr);
    exit (-1);                                             
  }

  /* read all the data */
  while (1) {
    /* read next line */
    fgets(buffer, BUFSIZE, DATA);
    if (feof(DATA)) break;

    /* stop if blank line */
    if (isspace(buffer[0]))
      break;

    /* read x and expected values for f(x) and f'(x) */
    ptr_buffer = buffer;
    x = strtod( ptr_buffer, &ptr_buffer );
    fx_exp = strtod( ptr_buffer, &ptr_buffer );
    dfx_exp = strtod( ptr_buffer, &ptr_buffer );

    /* compute function values for function tree */
    fx_obs = unur_distr_cont_eval_pdf(x,distr);
    dfx_obs = unur_distr_cont_eval_dpdf(x,distr);

    /* compare */




    fprintf(TESTLOG,"%g %g %g\n", x,fx_obs,dfx_obs);

  }

  /* free memory */
  unur_distr_free(distr);

  return 1;
#undef BUFSIZE
} /* end of _unur_test_function() */

/*---------------------------------------------------------------------------*/
