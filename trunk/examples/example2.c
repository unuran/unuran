/* my second UNURAN program example2.c */

#include <unuran.h>

/* desired cdf  -- piecewise parabolic */
double mycdf(double x, UNUR_DISTR *dummy)
{
   if (x <= 0.0)
     return ( 0.0 );
   else if (x <= 1.0)
     return ( x*x/3.0 );
   else if (x <= 3.0)
     return( ((-x+6.0)*x-3.0)/6.0 );
   else
     return ( 1.0 );
}

int main()
{
  double x;
  int    i;

  UNUR_DISTR *distr;  /* distribution */
  UNUR_PAR   *par;    /* parameter    */
  UNUR_GEN   *gen;    /* generator    */

  /* make distribution object and assign desired cdf (trig defined above) */
  distr = unur_distr_cont_new();
  unur_distr_cont_set_cdf(distr, mycdf);

  /* choose method and set parameters (Numerical INVersion) */
  par = unur_ninv_new(distr);

  /* make generator object */
  gen = unur_init(par);

  /* print 100 random numbers */
  for (i=0;i<100;i++) {
  x = unur_sample_cont(gen);
  printf("%f\n",x);
  }

  /* destroy distribution- and generator object */
  unur_distr_free(distr);
  unur_free(gen);

  exit (0);
}
