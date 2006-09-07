/* An example application testing the try_source.h interface. */
/* compile with : gcc -lm -ansi -pedantic try.c -o try */
/* on alpha machines add the compiler symbol -DALPHA */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "try_source.h"
 
int main() { 
  double x;

  x=1/0.; x+=0.; printf("%g\n",x);
  x=0/0.; x+=0.; printf("%g\n",x);

  x=sqrt(-1) ; printf("%g\n",x); 

  TRY { x=1/1; printf("1/1=%f\n",x); }
  CATCH(x) fprintf(stderr, "catching 1/1\n");

 x=17.; 
  TRY { x=0/0.; printf("0/0=%f\n",x); }
  CATCH(!(x==x)) x=0 ;

  printf(">>> %g\n",x);
  
 x=17.; 

  TRY { x=1/0.; printf("1/0=%f\n",x); }
  CATCH(!(x==x)) x=0;

  printf(">>> %g\n",x);

  TRY { x=sqrt(-1); printf("sqrt(-1)=%f\n",x); }
  CATCH(x>0) fprintf(stderr, "catching sqrt(-1)\n"); 

  return EXIT_SUCCESS;
}
