/* my second UNURAN program example2.c                             */

/* This method uses inversion for sampling. Therefore only
   a routine computing the cumulated distribution function
   is needed. Other methods might need the distribution function,
   the mode, ... 
   The section "Distributions" shows how to set them.              */   

#include <unuran.h>

/* --------------------------------------------------------------- */

/* CDF  -- piecewise parabolic                                     */
double mycdf(double x, UNUR_DISTR *dummy)
{
   if ( x <= 0.0 )
      return ( 0.0 );
   else if ( x <= 1.0 )
      return ( x*x/3.0 );
   else if ( x <= 3.0 )
      return( ((-x+6.0)*x-3.0)/6.0 );
   else
      return ( 1.0 );
}

/* --------------------------------------------------------------- */

int main()
{
  int    i;

  UNUR_DISTR *distr;  /* distribution                              */
  UNUR_PAR   *par;    /* parameter                                 */
  UNUR_GEN   *gen;    /* generator                                 */

  /* create a new distribution object (continuous distribution)
     this object is empty, to be of use, at least some values
     must be set (depending on the method used for sampling).      */
  distr = unur_distr_cont_new();

  /* Assign the cdf (defined above) to the distribution object     */
  unur_distr_cont_set_cdf(distr, mycdf);

  /* choose method and set parameters: Numerical INVersion         */
  par = unur_ninv_new(distr);

  /* make generator object                                         */
  gen = unur_init(par);
  
  if ( gen == NULL ){
     fprintf(stderr, "Error creating generation object\n");
     return (1);
  }

  /* print 100 random numbers                                      */
  for (i=0;i<100;i++)
     printf( "%f\n", unur_sample_cont(gen) );

  /* destroy distribution- and generator object                    */
  unur_distr_free(distr);
  unur_free(gen);

  return (0);
}
