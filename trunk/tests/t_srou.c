/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Tests                                                                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "t_unuran.h"

/*---------------------------------------------------------------------------*/

static const char method_name[] = "srou";

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

static int test_ok = TRUE;        /* all tests ok (boolean)                  */
static FILE *TESTLOG = NULL;      /* test log file                           */
static FILE *UNURANLOG = NULL;    /* unuran log file                         */

/*---------------------------------------------------------------------------*/

int main()
{ 
  char filename[128];

  /* open test log file */
  sprintf(filename,"%s_test.log",method_name);
  TESTLOG = fopen(filename,"w");
  if (TESTLOG == NULL) exit (-1);

  /* set output stream for unuran messages */
  sprintf(filename,"%s_unuran.log",method_name);
  UNURANLOG = fopen(filename,"w");
  if (UNURANLOG == NULL) exit (-1);
  unur_set_stream( UNURANLOG );
  
  /* close log files and exit */
  fclose(TESTLOG);
  fclose(UNURANLOG);
  exit( (test_ok) ? 0 : -1 );
}

/*---------------------------------------------------------------------------*/

