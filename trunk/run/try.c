/* An example application testing the try_source.h interface. */
/* compile with : gcc -lm -ansi -pedantic try.c -o try */
/* on alpha machines add the compiler symbol -DALPHA */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "try_source.h"
 
int main() { 
  double x;
  
  TRY { x=1/1; printf("1/1=%f\n",x); }
  CATCH(x) fprintf(stderr, "exception in 1/1\n");

  TRY { x=0/0; printf("0/0=%f\n",x); }
  CATCH(x) fprintf(stderr, "exception in 0/0\n");

  TRY { x=1/0; printf("1/0=%f\n",x); }
  CATCH(x) fprintf(stderr, "exception in 1/0\n");

  TRY { x=sqrt(-1); printf("sqrt(-1)=%f\n",x); }
  CATCH(x) {
    fprintf(stderr, "exception in sqrt(-1)\n");

    TRY { 
      x=1/1; printf("1/1=%f\n",x); 
      TRY { x=0/0; printf("0/0=%f\n",x); }
      CATCH(x) fprintf(stderr, "exception in 0/0\n");
    }
    CATCH(x) {
      fprintf(stderr, "exception in 1/1\n");
    }
  }

  return EXIT_SUCCESS;
}
