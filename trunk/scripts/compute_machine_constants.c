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

/** TODO:
#define MAXLGM 2.556348e305
**/

/*---------------------------------------------------------------------------*/

int main()
{

  printf("/** File automatically created by scripts/compute_machine_constants         **/\n\n");

  printf("/*****************************************************************************\n");
  printf(" *                                                                           *\n");
  printf(" *          unuran -- Universal Non-Uniform Random number generator          *\n");
  printf(" *                                                                           *\n");
  printf(" *****************************************************************************/\n");

  printf("\n");

  printf("/*---------------------------------------------------------------------------*/\n");
  printf("#ifndef __SOURCE_FP_CONST_H_SEEN\n");
  printf("#define __SOURCE_FP_CONST_H_SEEN\n");
  printf("/*---------------------------------------------------------------------------*/\n");

  printf("/* the machine roundoff error */\n");  
  printf("#define MACHEP  %.30g\n\n",DBL_EPSILON/2.);

  printf("/* largest argument for exp() */\n");
  printf("#define MAXLOG  %.30g\n\n",log(DBL_MAX)); 

  printf("/* smallest argument for exp() without underflow */\n");
  printf("#define MINLOG  %.30g\n\n",log(DBL_MIN)); 

  printf("/* the maximal number that pow(x,x-0.5) has no overflow */\n");
  printf("/* we use a (very) conservative portable bound          */\n");
  printf("#define MAXSTIR  %.30g\n\n", log(DBL_MAX) / log(log(DBL_MAX)) );

  printf("/*---------------------------------------------------------------------------*/\n");
  printf("#endif  /* __SOURCE_FP_CONST_H_SEEN */\n");
  printf("/*---------------------------------------------------------------------------*/\n");

  exit (0);
}

/*---------------------------------------------------------------------------*/
