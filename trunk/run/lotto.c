/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Lotto number generator : draw P distinct numbers from {1, 2, ..., N}     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#include <unuran.h>
#include <unur_struct.h>
#include <unur_source.h>

int N=45;
int P=6;
int SEED=163;

/*---------------------------------------------------------------------------*/

void printarray(int *a) {
  int  p;
  
  for (p=0; p<P; p++) {
    printf("  % d", a[p]);
  }
  printf("\n");

}


/*---------------------------------------------------------------------------*/

int main()
{

  UNUR_GEN *gen=NULL;
  int *a=NULL;
  char api_str[80];
  
  int i,p ;
  int flag_equal;
  double x;
  
  time_t now;
  struct tm t;

 
  unur_set_default_debug(UNUR_DEBUG_OFF);

  a=malloc(P*sizeof(int));

  time(&now);
  t = *localtime(&now);

  SEED = (int) (11723 +  SEED * (t.tm_sec*t.tm_min+(t.tm_sec*2.1)*t.tm_mon ));
  sprintf(api_str, "uniform(0,%d) & urng=mt19937(%d)", N, SEED); 
  gen=unur_str2gen(api_str);

  if (gen==NULL) return 0; 
  
  
  for (i=0; i<7+4*t.tm_sec; i++) {
    x=unur_sample_cont(gen);
  }


  p=0;
  while (p<P) {
    x=unur_sample_cont(gen);
    a[p] = 1 + (int) (x);

    flag_equal=0;
    for (i=0; i<p; i++) {
      if (a[i]==a[p]) flag_equal=1; 
    }
    if (flag_equal==0) p++;
  }

  printarray(a);
  
  unur_free(gen);
  
  if (a) free(a); 
    
  return 0;

}


