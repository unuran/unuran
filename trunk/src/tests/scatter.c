
/*  #error method UNIF has changed! */


/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      scatter.c                                                    *
 *                                                                           *
 *   make a scatter plot                                                     *
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
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

#define PLOT_DOT_SIZE   "0.01"
#define SCATTER_UNIFORM "uniform.scatterplot"

/*---------------------------------------------------------------------------*/
static char test_name[] = "Scatter";
/*---------------------------------------------------------------------------*/

static int _unur_make_uniform_scatter( int start, int skip );

/*---------------------------------------------------------------------------*/
/* baby generator                                                            */

#define UNUR_BABYGEN_PERIOD  1024 

#if UNUR_URNG_TYPE == UNUR_URNG_POINTER 
static double _unur_babygen( void );
#endif
static UNUR_URNG *_unur_get_babygen( void );

/*---------------------------------------------------------------------------*/

int
unur_make_scatterplot( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make scatterplot of generated numbers                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator                                    */
     /*   cdf    ... pointer to c.d.f. of distribution                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if successful                                                */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#if ( ( (UNUR_URNG_TYPE == UNUR_URNG_POINTER) || \
        (UNUR_URNG_TYPE == UNUR_URNG_PRNG) )  && \
        defined( UNUR_ENABLE_LOGGING ) )
/*---------------------------------------------------------------------------*/
{
#define DISTR   gen->distr.data.cont

  static int can_run_plotting_program = 1;  /* store failure */

  double Fl, Fr, Fdelta;  /* value of cdf (at left and right boundary point) */
  UNUR_FUNCT_CONT *cdf;                     /* pointer to c.d.f. */
  char *scatter_filename;                   /* name of scatter files */
  FILE *scatter;                            /* file handle for scatter files */
  char *call_graph;                         /* string for system call        */
  int len_scatter_filename, len_call_graph; /* string lengths                */
  double store, new;                        /* store generated points        */
  int n_urn;                            /* count number of generated points  */
  int start, skip;                      /* parameter for subsequence of urng */
  int error;                            /* exit code of system call          */  

  UNUR_URNG *urng_bak, *urng_aux_bak;   /* for saving URNG of generator      */
  UNUR_URNG *urng_babygen;         /* pointer to uniform baby RNG            */

  /* check arguments */
  _unur_check_NULL(test_name,gen,0);
  /* we do not check magic cookies here */

  /* c.d.f. required */
  cdf = DISTR.cdf;
  if (cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"c.d.f. required");
    return 0;
  }

  if (!_unur_gen_is_cont(gen)) {
    /* cannot make scatter plot */
    _unur_error(test_name,UNUR_ERR_GENERIC,"Not implemented for this type of generator");
    return 0;
  }

  /* has system call for plotting program been failure */
  if (!can_run_plotting_program)
    return 0;

  /* compute Fl and Fr */
  Fl = (DISTR.domain[0] <= -INFINITY) ? 0. : cdf(DISTR.domain[0],&(gen->distr));
  Fr = (DISTR.domain[1] >=  INFINITY) ? 1. : cdf(DISTR.domain[1],&(gen->distr));
  Fdelta = Fr - Fl;

  /* Fr - Fl <= 0. is a fatal error */
  if (Fdelta <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
    return -1;
  }

  /* which subsequence of the underlying urng should be used 
     for comparation? */
  switch (gen->method) {
  case UNUR_METH_AROU:
  case UNUR_METH_SROU:
  case UNUR_METH_TABL:
    start = 1;
    skip = 1;
    break;
  case UNUR_METH_SSR:
  case UNUR_METH_TDR:
  case UNUR_METH_UTDR:
    start = 1;
    skip = 2;
    break; 
  default: /* unknown ! */
    _unur_error(test_name,UNUR_ERR_GENERIC,"method unknown!");
    return 0;
  }

  /* make scatterplot for underlying uniform generator */
  if (! _unur_make_uniform_scatter(start,skip) ) {
    can_run_plotting_program = 0;
    return 0;
  }
 
  /* save generators */
  urng_bak = gen->urng;
  urng_aux_bak = gen->urng_aux;

  /* store uniform random number generator */
  urng_babygen = _unur_get_babygen();
  if (urng_babygen == NULL) return 0;
  unur_chg_urng(gen, urng_babygen);

  /* we do not change the auxilliary URNG */
  unur_chg_urng_aux(gen, urng_aux_bak);

  /* make string for scatter file name */
  len_scatter_filename = strlen(gen->genid) + strlen(".scatterplot") + 2;
  scatter_filename = _unur_malloc( len_scatter_filename * sizeof(char) );
  sprintf(scatter_filename,"%s.scatterplot",gen->genid);

  /* open file for scatter plot */
  scatter = fopen(scatter_filename,"w");
  if (scatter == NULL) return 0;

  /* set starting value */
  store = cdf( unur_sample_cont(gen), &(gen->distr));

  /* sampling */
  for (n_urn = 0; n_urn < UNUR_BABYGEN_PERIOD; ++n_urn ) {

    /* invoke generator and scale sample according to cdf */
    new = cdf( unur_sample_cont(gen), &(gen->distr));

    /* write pair samples into scatter file */  
    fprintf(scatter,"%f %f\n",store, new); 
    
    /* move frame for pair one position ahead */     
    store = new;
  }

  /* close scatter file */
  fclose(scatter);
  
  /* restore uniform random number generator */
  unur_chg_urng(gen, urng_bak);

  /* make plots */
  len_call_graph = ( strlen("graph -T X -C -m0 -S 16 %s %s -L%s -C -m-2 -S 16 %s %s")
		     + 2*strlen(PLOT_DOT_SIZE) +100
		     + strlen(SCATTER_UNIFORM) + strlen(scatter_filename)
		     + strlen(gen->genid) );
  call_graph = _unur_malloc( len_call_graph * sizeof(char) );
/*    sprintf(call_graph,"graph -T X -C -m-4 -S 16 %s %s -L%s -C -m-1 -S 16 %s %s ", */
/*  	  PLOT_DOT_SIZE,  */
/*  	  SCATTER_UNIFORM, */
/*  	  gen->genid,  */
/*  	  PLOT_DOT_SIZE,  */
/*  	  scatter_filename);  */
  sprintf(call_graph,"graph -T X -C -m-4 -L%s -C -m-1 -S 16 %s %s ",
	  gen->genid, 
	  PLOT_DOT_SIZE, 
	  scatter_filename); 

  /* run system call */
  error = system(call_graph);
  if (error) {
    _unur_warning(test_name,UNUR_ERR_GENERIC,"Cannot run \"graph\"");
    can_run_plotting_program = 0;
  }

  /* free memory */
  free(scatter_filename);
  free(call_graph);
    
  /* o.k. */
  return 1;

#undef DISTR
} /* end of unur_make_scatterplot() */
/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
{
  _unur_error(test_name,UNUR_ERR_GENERIC,"Cannot make scatter plot.\n Recompile with different UNUR_URNG_TYPE!\n Set flag UNUR_ENABLE_LOGGING");
  return -1;
} /* end of unur_make_scatterplot() */
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

static int
_unur_make_uniform_scatter( int start, int skip )
     /*----------------------------------------------------------------------*/
     /* make scatterplot for uniform random number  generator                */
     /*                                                                      */
     /* parameter:                                                           */
     /*   start   ... starting point of sequence of baby generator           */
     /*   skip    ... skip for subsequence                                   */
     /*                                                                      */
     /*   return:                                                            */
     /*     1 ... on success                                                 */
     /*     0 ... on error                                                   */
     /*----------------------------------------------------------------------*/
{
  FILE *scatter;                         /* file handle to scatter_filename  */
  double store, new;                     /* store generated points           */
  int n_urn;                             /* count number of generated points */
  struct unur_gen *gen;                  /* pointer to generator object      */
  UNUR_URNG *urng_babygen;               /* pointer to uniform RNG           */

  /* make string for scatter file name */

  /* open file for scatter plot */
  scatter = fopen(SCATTER_UNIFORM,"w");
  if (scatter == NULL)  return 0;

  /* make generator object for uniform baby generator */
  gen = unur_init( unur_unif_new() );
  _unur_check_NULL(test_name,gen,0 );

  /* get pointer to baby generator */
  urng_babygen = _unur_get_babygen();
  if (urng_babygen == NULL) return 0;
  unur_chg_urng(gen,urng_babygen);

  /* get starting value */
  store = unur_sample_cont(gen);

  for (n_urn = 0; n_urn < UNUR_BABYGEN_PERIOD; ++n_urn ) {

    /* invoke generator */
    new = unur_sample_cont(gen);
    
    /* write sampled tuple into scatter file */
    fprintf(scatter,"%f %f\n",store, new); 

    /* move frame for tuple one position ahead */     
    store = new;
  }

  /* close the file */
  fclose(scatter);

  /* destroy generator object */
  _unur_free(gen);

  /* o.k. */
  return 1;
   
}/* end of _unur_make_uniform_scatter() */

/*--------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   baby generator                                                          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_POINTER 
/*---------------------------------------------------------------------------*/

static int u = 0;   /* state variable of baby generator                      */

/*---------------------------------------------------------------------------*/
double
_unur_babygen( void )
     /*----------------------------------------------------------------------*/
     /* baby generator (LCG with extreme short period)                       */
     /*                                                                      */
     /* return:                                                              */
     /*   uniform random number                                              */
     /*----------------------------------------------------------------------*/
{
  u = (869*u+1) % 1024;
  return u/1024.0;
} /* end of _unur_babygen() */

/*--------------------------------------------------------------------------*/

UNUR_URNG *
_unur_get_babygen( void )
     /*----------------------------------------------------------------------*/
     /* evoke baby generator and reset seed                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to uniform random number generator                         */
     /*----------------------------------------------------------------------*/
{
  static UNUR_URNG *urngen = NULL;

  /* set babygen */
  if (urngen == NULL)
    urngen = _unur_babygen;

  /* reseed the generator */
  u = 0;

  return urngen;
} /* end of _unur_get_babygen */

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

UNUR_URNG *
_unur_get_babygen( void )
     /*----------------------------------------------------------------------*/
     /* evoke baby generator and reset seed                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to uniform random number generator                         */
     /*----------------------------------------------------------------------*/
{
  static UNUR_URNG *urngen = NULL;

  /* set babygen */
  if (urngen == NULL) {
    urngen = prng_new("LCG(1024,869,1,0)");
    if( urngen == NULL ) {
      /* some parameters invalid! */
      _unur_error("prng",UNUR_ERR_NULL,"Cannot set baby generator");
      return NULL;
    }
  }

  /* reseed the generator */
  prng_seed(urngen,1);

  return urngen;
} /* end of _unur_get_babygen */

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/

static UNUR_URNG_TYPE
_unur_get_babygen( void )
     /*----------------------------------------------------------------------*/
     /* not implemented for this choice of UNUR_URNG_TYPE !!               */
     /*----------------------------------------------------------------------*/
{
  return NULL;
} /* end of _unur_get_babygen */

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/



