/* Test for signalling FP exceptions */

#define ALPHA 0 /* set to 1 on alpha machine */

#include <stdio.h>
#include <signal.h>
#include <setjmp.h>

#if !ALPHA
#include <fenv.h> /* not available on alphas */
#endif 
 
jmp_buf jump_mark; /* Address for long jump to jump to */

/* our floating point exception handler routine */
void fphandler( int sig )
{
   fprintf(stderr, "invoking fphandler ...\n");
   
   /* Restore calling environment and jump back to setjmp. */
   /* Return -1 so that setjmp will return false for conditional test. */
   longjmp( jump_mark, -1 );
}

/*---------------------------------------------------------------------*/

/* f(x,y) = x/y */
double f(double x, double y) {

  double fx;
  
  int jmpret;
  
#if !ALPHA  
  /* not available on alphas */ 
  feenableexcept(FE_DIVBYZERO); 
  feenableexcept(FE_ALL_EXCEPT); 
#endif

  fprintf(stderr, "----------------------------------\n");
  fprintf(stderr, "invoking f1(%d, %d)\n", (int) x, (int) y);
  
#if !ALPHA
  /* not available on alphas */ 
  fprintf(stderr, "exceptions enabled: %x\n", fegetexcept()); 
#endif  
  
  jmpret = setjmp( jump_mark );

  if( jmpret == 0 ) {
    fx=x/y ; fprintf(stderr, "%d/%d = %f\n", (int) x, (int) y, fx);
  }
  else { 
    fx=0;
    fprintf(stderr, 
            "%d/%d ... could not be calculated (jmpret=%d)\n", 
            (int) x, (int) y, jmpret);  
  }
  
  return fx;
}

/*---------------------------------------------------------------------*/

int main()
{
  if( signal( SIGFPE, fphandler ) == SIG_ERR ) {
    fprintf( stderr, "Couldn't set SIGFPE\n" );
  }

  printf("f(%d,%d)=%f\n", 0,0,f(0,0)); 
  printf("f(%d,%d)=%f\n", 0,1,f(0,1)); 
  printf("f(%d,%d)=%f\n", 1,0,f(1,0)); 
  printf("f(%d,%d)=%f\n", 1,1,f(1,1)); 

  fprintf(stderr, "OK\n");

  return 0;
}

/*---------------------------------------------------------------------*/
