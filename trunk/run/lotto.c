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

void printarray(char *info, int *a) {
  int  p;
  
  printf("%s", info);

  for (p=0; p<P; p++) {

    printf("  % d", a[p]);
  }
  
  printf("\n");

}


int main(int argc, char *argv[])
{

  UNUR_GEN *gen=NULL;

  int *a;
  char api_str[40];
  char info[20];
  
  int i,p ;
  int flag_equal;
  double x;
  
  char c;
  
  time_t now;
  struct tm t;

  /* read options */
  while ((c = getopt(argc, argv, "p:n:h")) != -1) {
    switch (c) {
    case 'p':     
      P=atol(optarg);
      break;
    case 'n':     
      N=atol(optarg);
      break;
    
    case 'h':     
      printf("options\n" );
      printf(" -n N         : default = %d  \n", N );
      printf(" -p P         : default = %d  \n", P );
      
      return 0;
      break;
    default:
      break;
    }
  }
  
  if (N>1 && N<100 && 
      P>1 && P<N) {
        a=malloc(P*sizeof(int));
  }    
  else {
    return 0;
  }

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

  printarray("", a);
  
  unur_free(gen);
  
  if (a) free(a); 
    
  return 0;

}

/*---------------------------------------------------------------------------*/

