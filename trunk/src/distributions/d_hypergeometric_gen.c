/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_hypergeometric_gen.c                                       *
 *                                                                           *
 *   Special generators for Hypergeometric distribution                      *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.dstd        /* data for parameter object         */
#define GEN       gen->data.dstd        /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params   8     /** TODO:  maximal number of parameters for generator */
#define MAX_gen_iparams  3     /** TODO:  maximal number of parameters for generator */

/* parameters */
#define N  (DISTR.params[0])
#define M  (DISTR.params[1])
#define n  (DISTR.params[2])

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

#ifdef HAVE_UNUR_SF_LN_FACTORIAL
inline static int hypergeometric_hruec_init( struct unur_gen *gen );
int _unur_stdgen_sample_hypergeometric_hruec( struct unur_gen *gen );
#endif

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_hypergeometric_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Hypergeometric distribution         */
     /* if gen == NULL then only check existance of variant.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* HRUEC  method */
#ifdef HAVE_UNUR_SF_LN_FACTORIAL
     _unur_dstd_set_sampling_routine( par,gen,_unur_stdgen_sample_hypergeometric_hruec );
     return hypergeometric_hruec_init( gen );
#else
     return 0;
#endif

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_hypergeometric_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

#ifdef HAVE_UNUR_SF_LN_FACTORIAL

#define flogfak(k) _unur_sf_ln_factorial(k)
/*---------------------------------------------------------------------------*/

/******************************************************************
 *                                                                *
 *   Hypergeometric Distribution - Ratio of Uniforms/Inversion    *
 *                                                                *
 ******************************************************************
 *                                                                *
 * Ratio of Uniforms combined with Inversion for                  *
 * sampling from Hypergeometric distributions with parameters     *
 * N, M and n. The algorithm is valid for M <= N/2, n <= N/2.     *
 * Otherwise parameters (at the beginning of the algorithm) and   *
 * random numbers k are adapted in function h_util().             *
 * For mode m < 5 Inversion is applied:                           *
 * The random numbers are generated via sequential search,        *
 * starting at the lowest index k=0. The cumulative probabilities *
 * are avoided by using the technique of chop-down.               *
 * For mode  m >=5  Ratio of Uniforms is employed:                *
 * A table mountain hat function h(x) with optimal scale          *
 * parameter s for fixed location parameter  a = mu+1/2  is used. *
 * If the mode m <= 20 and the candidate k is near the mode       *
 * f(k) is computed recursively starting at the mode  m.          *
 *                                                                *
 ******************************************************************
 *                                                                *
 * FUNCTION:    - hruec samples a random number from the          *
 *                Hypergeometric distribution with parameters     *
 *                N (number of red and black balls), M (number    *
 *                of red balls) and n (number of trials)          *
 *                valid for N >= 2, M,n <= N.                     *
 * REFERENCE:   - E. Stadlober (1989): Sampling from Poisson,     *
 *                binomial and hypergeometric distributions:      *
 *                ratio of uniforms as a simple and fast          *
 *                alternative, Bericht 303, Math. Stat. Sektion,  *
 *                Forschungsgesellschaft Joanneum, Graz.          *
 * SUBPROGRAMS: - flogfak(k)  ... log(k!) with integer k     *
 *                                                                *
 * Implemented by R.Kremer 1990, revised by P.Busswald, July 1992 *
 ******************************************************************/

#define delta(k) (flogfak(k)+flogfak(Mc-k)+flogfak(nc-k)+flogfak(NMn+k))



#define b       (GEN.gen_iparam[0])
#define m       (GEN.gen_iparam[1])
#define NMn     (GEN.gen_iparam[2])

#define NMnp    (GEN.gen_param[0])
#define Np      (GEN.gen_param[1])
#define Mp      (GEN.gen_param[2])
#define np      (GEN.gen_param[3])
#define g       (GEN.gen_param[4])
#define a       (GEN.gen_param[5])
#define h       (GEN.gen_param[6])
#define p0      (GEN.gen_param[7])


#define ln2 0.69314718055994531


/*---------------------------------------------------------------------------*/

inline static int
hypergeometric_hruec_init( struct unur_gen *gen )
{
 int N_half,Mc,nc,k,k1,i,bh;
 double x,u,f,lf,p,q,c,my;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DSTD_GEN,0);

  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = MAX_gen_params;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
    GEN.n_gen_iparam = MAX_gen_iparams;
    GEN.gen_iparam = _unur_malloc(GEN.n_gen_param * sizeof(int));
  }

  /* -X- setup code -X- */

  N_half=N/2.0;                      /* Preparations of the parameters */
 if (M<=N_half) Mc=M;else Mc=N-M;   /* if M<=N/2 M is replaced by N-M */
 if (n<=N_half) nc=n;else nc=N-n;   /* if n<=n/2 n is replaced by N-n */
 	                                     /* Set-up */
	 Np=(double)N;
	 Mp=(double)Mc;
	 np=(double)nc;
	 NMn=N-Mc-nc;
	 NMnp=Np-Mp-np;
	 p=Mp/Np;
	 q=1.0-p;
	 my=np*p;
	 bh=min(nc,Mc);
	 m=(int) ((np+1.0)*(Mp+1.0)/(Np+2.0));       /* mode */
	 if (m<5)
	{
	 c=my+10.0*sqrt(my*q*(1.0-np/Np));     /* Set-up for Inversion */
	 b=min(bh,(int)c);                /* safety-bound */
	 p0=exp(flogfak(N-Mc)+flogfak(N-nc)-flogfak(NMn)-flogfak(N));
	}
	 else
	{
	 a=my+0.5;                      /* Set-up for Ratio of Uniforms */
	 c = sqrt(2.0*a*q*(1.0-np/Np));
	 b=min(bh,(int)(a+7.0*c));       /* safety-bound */
	 g=delta(m);
	 k1=(int)(a-c);
	 x=(a-k1-1.0)/(a-k1);
	 if((np-k1)*(p-(double)k1/Np)*x*x > (k1+1)*(q-(np-k1-1.0)/Np)) k1++;
	 h=(a-k1)*exp(0.5*(g-delta(k1))+ln2);    /* h=2*s */
	}
	 


  /* -X- end of setup code -X- */

  return 1;

} /* end of hypergeometric_hruec_init() */

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_hypergeometric_hruec( struct unur_gen *gen )
{
  /* -X- generator code -X- */
 int N_half,Mc,nc,k,k1,i,bh, h_util(int iN,int iM,int in,int ik);
 double x,u,f,lf,p,q,c,my;
  

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DSTD_GEN,0);


 if (m<5)                                     /* Inversion/Chop-down */
	 {
	  double pk;

	  k=0;
	  pk=p0;
	  u=uniform();
	  while (u>pk)
		{
		 ++k;
		 if (k>b)
			 {
		u=uniform();
		k=0;
		pk=p0;
			 }
		 else
			 {
		u-=pk;
		pk*=((Mp-k+1.0)*(np-k+1.0))/((double)k*(NMnp+k));
			 }
		}
	  return (h_util(N,M,n,k));
	 }
 for (;;)                                     /* Ratio of Uniforms */
	 {
	  do
		 {
	u=uniform();
	x=a+h*(uniform()-0.5)/u;
		 }
	  while (x < 0 || ((k=(int)x) > b)); /* check, if k is valid candidate */
	  if (m<=20 || labs(m-k)<=15)
		 {                                     /* compute f(k) recursively */
	f=1.0;
	if (m<k)
	  {
		for (i=m+1;i<=k;i++) f*=((Mp-i+1.0)*(np-i+1.0))/((double)i*(NMnp+i));
		if (u*u <= f) break;              /* f - f(k), u^2<=f */
	  }
	else
	  {
		for (i=k+1;i<=m;i++) f*=((Mp-i+1.0)*(np-i+1.0))/((double)i*(NMnp+i));
		if (u*u*f <= 1.0) break;          /* f - 1/f(k), u^2<=f */
	  }
	}
	  else
		 {
	lf=g-delta(k);                       /* lf - ln(f(k)) */
	if ( u * (4.0 - u) - 3.0 <= lf) break;     /* lower squeeze */
	if (u*(u-lf) <= 1.0)                       /* upper squeeze */
		if (2.0*log(u) <= lf) break;         /* final acceptance */
		 }
	}
 return (h_util(N,M,n,k));
} /* end of _unur_stdgen_sample_hypergeometric_hruec() */


#undef b    
#undef m   
#undef NMn 

#undef NMnp
#undef Np  
#undef Mp  
#undef np  
#undef g   
#undef a   
#undef h   
#undef p0  


#undef ln2 0.69314718055994531


#undef N
#undef M
#undef n
#undef delta


int h_util(int N,int M,int n,int k)
         /* Transformation to variate k from H(N,M,n) */

{
 int N_half;

 N_half=N/2.0;
 if (n<=N_half)
	 {
	  if (M<=N_half) return(k);   /* no replacements */
	  else return(n-k);           /* M has been replaced by N-M, therefore */
	 }                            /* k has to be replaced by n-k           */
 else
	 {
	  if (M<=N_half) return(M-k); /* n h.b.r. by N-n, therefore k by M-k */
	  else return(M-N+n+k);       /* M h.b.r. by N-M and n by N-n,       */
	 }                            /* therefore k by M-N+n+k              */
}

  

  /* -X- end of generator code -X- */



/*---------------------------------------------------------------------------*/
#endif   /*  HAVE_UNUR_SF_LN_FACTORIAL  */
/*---------------------------------------------------------------------------*/

















