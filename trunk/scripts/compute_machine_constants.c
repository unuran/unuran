/*****************************************************************************
 *                                                                           *
 *          unuran -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *  Compute some machine constants and print to stdout                       *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*---------------------------------------------------------------------------*/

int main()
{

#if 0
  printf("/* Include into file `src/unuran_config.h' */\n\n");

  printf("/* maximum number for that exp() can be computed without overflow */\n");
  printf("#define UNUR_MAX_ARG_EXP  %20.16f\n\n",log(DBL_MAX)); 

  printf("#define MACHEP  %.21g\n\n",DBL_EPSILON/2); 

#endif

  printf("log(DBL_MAX) = %.21g\n\n",log(DBL_MAX)); 
  printf("log(DBL_MIN) = %.21g\n\n",log(DBL_MIN)); 

  exit (0);
}

/*---------------------------------------------------------------------------*/
