/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_codegen.c                                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Make generator code for TDR.                                         *
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

#include <source_unuran.h>
#include "PDFgen_source.h"

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define TDR_VARMASK_T          0x000fu   /* indicates transformation         */
#define TDR_VAR_T_SQRT         0x0001u   /* T(x) = -1/sqrt(x)                */
#define TDR_VAR_T_LOG          0x0002u   /* T(x) = log(x)                    */
#define TDR_VAR_T_POW          0x0003u   /* T(x) = -x^c                      */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define GEN       gen->data.tdr         /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

int
_unur_tdr_ps_codegen( struct unur_gen *gen, FILE *out, const char *distr_name )
{
  const char hline[] = "\n/* ------------------------- */\n/* %-25s */\n/* ------------------------- */\n\n";

  struct unur_tdr_interval *iv;
  int i,j;

  /* names for function */
  char *pdf, *rand;

  /* check arguments */
  _unur_check_NULL("unurgen",gen, 0);
  COOKIE_CHECK(gen,CK_TDR_GEN,0);

  /* name of PDF function and sampling routine */
  if (distr_name == NULL) 
    distr_name = unur_distr_get_name( &(gen->distr) );
  pdf = _unur_malloc((5+strlen(distr_name)) * sizeof(char));
  sprintf(pdf,"pdf_%s",distr_name);
  rand = _unur_malloc((6+strlen(distr_name)) * sizeof(char));
  sprintf(rand,"rand_%s",distr_name);

  /* PDF */
  fprintf(out,hline,"PDF");
  if (! _unurgen_C_PDF(&(gen->distr),out,pdf)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"Cannot make PDF");
    free(pdf); free(rand);
    return 0;
  }

  /* sampling routine */
  fprintf(out,hline,"Sampling");
  fprintf(out,"double %s(void)\n{\n",rand);

  /* constants */
  fprintf(out,"\t/* data */\n");

  /* guide table */
  fprintf(out,"\tconst int guide_size = %d;\n", GEN.guide_size);
  fprintf(out,"\tconst int guide[%d] = {\n\t\t", GEN.guide_size);
  iv = GEN.iv;
  j = 0;
  for (i=0; i<GEN.guide_size; i++) {
    while (GEN.guide[i] != iv) {
      iv=iv->next;
      j++;
    }
    if (iv == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 0;
    }
    if (i==0)
      fprintf(out,"%d",j);
    else
      fprintf(out,(i%16)?", %d":",\n\t\t%d",j);
  }      
  fprintf(out," };\n");

  /* hat function */
  fprintf(out,"\tconst double Atotal = %.20g;\n", GEN.Atotal);

  fprintf(out,"\tconst struct {\n");
  fprintf(out,"\t\tdouble x;\n");
  fprintf(out,"\t\tdouble fx;\n");
  fprintf(out,"\t\tdouble Tfx;\n");
  fprintf(out,"\t\tdouble dTfx;\n");
  fprintf(out,"\t\tdouble sq;\n");
  fprintf(out,"\t\tdouble Acum;\n");
  fprintf(out,"\t\tdouble Ahatr;\n");
  fprintf(out,"\t} iv[%d] = {\n",GEN.n_ivs);
  
  for (iv=GEN.iv; iv->next!=NULL; iv=iv->next) {
    fprintf(out,(iv==GEN.iv)?"\t\t{":",\n\t\t{");
    fprintf(out," %.20g",iv->x);
    fprintf(out,", %.20g",iv->fx);
    fprintf(out,", %.20g",iv->Tfx);
    fprintf(out,", %.20g",iv->dTfx);
    fprintf(out,", %.20g",iv->sq);
    fprintf(out,", %.20g",iv->Acum);
    fprintf(out,", %.20g",iv->Ahatr);
    fprintf(out," }");
  }
  fprintf(out," };\n");

  fprintf(out,"\n");

  fprintf(out,"\t/* code */\n");

  /* declare variables */
  fprintf(out,"\tint I;\n");          /* random interval */
  fprintf(out,"\tdouble U;\n");       /* uniform random number */
  fprintf(out,"\tdouble V;\n");       /* random number for rejection */
  fprintf(out,"\tdouble X;\n");       /* non-uniform random variate */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\tdouble t;\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"\tdouble Thx;\n");
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  } /* end switch */
  fprintf(out,"\n");

  /* main loop */
  fprintf(out,"\twhile (1) {\n");

  /* sample from U(0,1) */
  fprintf(out,"\t\tU = uniform();\n");

  /* look up in guide table and search for interval */
  fprintf(out,"\t\tI =  guide[(int) (U * guide_size)];\n");
  fprintf(out,"\t\tU *= Atotal;\n");
  fprintf(out,"\t\twhile (iv[I].Acum < U) I++;\n");

  /* reuse of uniform random number */
  fprintf(out,"\t\tU -= iv[I].Acum - iv[I].Ahatr;\n");
  /* result: U in (-A_hatl, A_hatr) */

  /* generate from hat distribution */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\t\tt = iv[I].dTfx * U / iv[I].fx;\n");
    fprintf(out,"\t\tif (fabs(t) > 1.e-8)\n");
    fprintf(out,"\t\t\tX = iv[I].x + log(t + 1.) * U / (iv[I].fx * t);\n");
    fprintf(out,"\t\telse\n");
    fprintf(out,"\t\t\tX = iv[I].x + U / iv[I].fx * (1 - t/2.);\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"\t\tX = iv[I].x + U / iv[I].fx;\n");
    fprintf(out,"\t\tif (iv[I].dTfx != 0.)\n");
    fprintf(out,"\t\t\tX /= (1. - iv[I].Tfx * iv[I].dTfx * U);\n");
    break;
  } /* end switch */

  /* accept or reject */
  fprintf(out,"\t\tV = uniform();\n");

  /* squeeze acceptance */
  fprintf(out,"\t\tif (V <= iv[I].sq) return X;\n");

  /* evaluate hat at X:
     get uniform random number between 0 and hat(X) */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\t\tV *= iv[I].fx * exp(iv[I].dTfx * (X - iv[I].x));\n");
    break;
  case TDR_VAR_T_SQRT:
     /* transformed hat at X */ 
    fprintf(out,"\t\tThx = iv[I].Tfx + iv[I].dTfx * (X - iv[I].x);\n");
    fprintf(out,"\t\tV /= Thx*Thx;\n");
    break;
  } /* end switch */

  /* main rejection */
  fprintf(out,"\t\tif (V <= %s(X)) return X;\n",pdf);

  /* end of loop */
  fprintf(out,"\t}\n");

  /* end of function */
  fprintf(out,"}\n");

  fprintf(out,hline,"End");

  /* free memory */
  free(pdf); free(rand);

  /* o.k. */
  return 1;

} /* end of _unur_tdr_ps_codegen() */

/*****************************************************************************/

#if 0
struct unur_tdr_interval {

  double  x;                    /* (left) construction point (cp)            */
  double  fx;                   /* value of PDF at cp                        */ 
  double  dTfx;                 /* derivative of transformed PDF at cp       */
  double  sq;                   /* slope of transformed squeeze in interval  */
  double  Acum;                 /* cumulated area of intervals               */
  double  Ahatr;                /* area below hat on right side              */
};
#endif


/*****************************************************************************/
#if 0

double
_unur_tdr_ps_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (proportional squeeze)                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   log(x):                                                            */
     /*                                                                      */
     /*   hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                        */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                       */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  double U, V, X;
  double fx, Thx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U( Umin, Umax ) */
    U = GEN.Umin + _unur_call_urng(urng) * (GEN.Umax - GEN.Umin);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum - iv->Ahatr;    /* result: U in (-A_hatl, A_hatr) */

    /* generate from hat distribution */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      if (iv->dTfx == 0.)
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  /* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;

    case TDR_VAR_T_SQRT:
      if (iv->dTfx == 0.)
	X = iv->x + U /iv->fx;
      else {
	/* it would be less expensive to use:
	   X = iv->x + iv->Tfx/iv->dTfx * (1. - 1./(1. + iv->dTfx * iv->Tfx * U) )
	   however, this is unstable for small iv->dTfx */
	// 	X = iv->x + (iv->Tfx*iv->Tfx*U) / (1.-iv->Tfx*iv->dTfx*U);  
	X = iv->x + (U / iv->fx) / (1.-iv->Tfx*iv->dTfx*U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;

      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* accept or reject */
    V = _unur_call_urng(urng);

    /* squeeze rejection */
    if (V <= iv->sq)
      	return X;

    /* evaluate hat at X:
       get uniform random number between 0 and hat(X) */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      V *= iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      V *= 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
    } /* end switch */

    /* evaluate PDF at X */
    fx = PDF(X);

    /* main rejection */
    if (V <= fx)
      return X;

    /* else reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }

} /* end of _unur_tdr_ps_sample() */

/*****************************************************************************/

double
_unur_tdr_ia_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (immediate acceptance)                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   log(x):                                                            */
     /*                                                                      */
     /*   hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                        */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                       */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  int use_ia;
  double U, V, X;
  double fx, hx, Thx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U(0,1) */
    U = _unur_call_urng(urng);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat,0) */

    /* check for region of immediate acceptance */
    if (U >= - iv->sq * iv->Ahat) {
      /* region of immediate acceptance */
      U /= iv->sq;
      use_ia = 1;
    }
    else {
      /* rejection from region between hat and squeeze */
      U = (U + iv->sq * iv->Ahat) / (1. - iv->sq);
      use_ia = 0;
    }
    /* result: U in (-A_hat,0) */

    /* U in (-A_hatl, A_hatr) */
    U += iv->Ahatr;

    /* generate from hat distribution */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      if (iv->dTfx == 0.)
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  /* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;

    case TDR_VAR_T_SQRT:
      if (iv->dTfx == 0.)
	X = iv->x + U /iv->fx;
      else {
	U *= iv->Tfx; /* avoid one multiplication */
	X = iv->x + (iv->Tfx * U) / (1. - iv->dTfx * U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* immedate acceptance */
    if (use_ia)
      return X;

    /* evaluate hat at X */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      hx = 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      /** TODO **/
      return 1.;
    } /* end switch */

    /* from now on we use the auxilliary generator
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

    /* rejection from region between hat and (proportional) squeeze */
    V = _unur_call_urng(urng);

    /* get uniform random number between squeeze(X) and hat(X) */
    V = (iv->sq + (1 - iv->sq) * V) * hx;

    /* evaluate PDF at X */
    fx = PDF(X);

    /* main rejection */
    if (V <= fx)
      return X;

    /* else reject and try again */
  }

} /* end of _unur_tdr_ia_sample() */

/*****************************************************************************/
#endif
