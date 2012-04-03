/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ptx_init.c                                                  *
 *                                                                           *
 *   Routines for initialization and deletion of generator objects.          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_ptx_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params  pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  struct unur_par *pin;
  double lfm;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_PTX ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_PTX_PAR,NULL);

  /* GIN */
  pin = unur_pinv_new(PAR->din);
  unur_pinv_set_u_resolution(pin,1.e-15);
  /** FIXME: copy URNG **/


  /* create a new empty generator object */    
  gen = _unur_ptx_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;


  /* GIN */
  gin = _unur_init(pin);
  if (gin==NULL) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"init of GIN failed");
      _unur_ptx_free(gen); return NULL;
  }    
  GEN->Tmin = GIN->bleft;
  GEN->Tmax = GIN->bright;


  /* compute rescaling factor for PDF */
  /* (only used when logPDF is given) */
  if (DISTR.logpdf != NULL &&
      (gen->distr->set & UNUR_DISTR_SET_MODE) &&
      ! (gen->variant & PTX_VARIANT_CDF) &&
      !_unur_FP_less(DISTR.mode,DISTR.domain[0]) &&
      !_unur_FP_greater(DISTR.mode,DISTR.domain[1]) ) {
    lfm = (DISTR.logpdf)(DISTR.mode,gen->distr);
    /* rescaling results in more evaluations of the logPDF, */
    /* when the logPDF is approximately 0.                  */
    /* so we only rescale the logPDF when it is too small.  */
    if (lfm < -3.)
      GEN->logPDFconstant = (DISTR.logpdf)(DISTR.mode,gen->distr);
  }

  /* check parameters */
  if (_unur_ptx_check_par(gen) != UNUR_SUCCESS) {
    _unur_ptx_free(gen); return NULL;
  }


#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_ptx_debug_init_start(gen);
#endif

  /* 1. Preprocessing:                                     */
  /*   find interval for computing Newton interpolation */
  if (_unur_ptx_preprocessing(gen) != UNUR_SUCCESS) {
    /* preprocessing failed */
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_ptx_debug_init(gen,FALSE);
#endif
    _unur_ptx_free(gen); return NULL;
  }
  
  /* compute table for Newton interpolation */
  if (_unur_ptx_create_table(gen) != UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_ptx_debug_init(gen,FALSE);
#endif
    _unur_ptx_free(gen); return NULL;
  }

  /* we do not need the table with CDF values any more. */
  /* thus we free the allocated memory.                 */
  _unur_lobatto_free(&(GEN->aCDF));
  _unur_lobatto_free(&(GEN->aTRX));

  /* make guide table */
  /*   _unur_ptx_make_guide_table(gen); */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_ptx_debug_init(gen,TRUE);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_ptx_init() */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_ptx_reinit( struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* re-initialize (existing) generator.                                  *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int rcode; */

/*   /\* check parameters *\/ */
/*   if ( (rcode = _unur_ptx_check_par(gen)) != UNUR_SUCCESS) */
/*     return rcode; */

/*   return UNUR_SUCCESS; */
/* } /\* end of _unur_ptx_reinit() *\/ */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_ptx_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_PTX_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_ptx_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_PTX_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_ptx_getSAMPLE(gen);
  gen->destroy = _unur_ptx_free;
  gen->clone = _unur_ptx_clone;
  /* gen->reinit = _unur_ptx_reinit; */

  /* copy parameters into generator object */
  GEN->order = PAR->order;            /* order of polynomial                 */
  GEN->u_resolution = PAR->u_resolution; /* maximal error in u-direction     */
  GEN->bleft_par  = PAR->bleft;          /* border of computational domain   */
  GEN->bright_par = PAR->bright;
  GEN->sleft  = PAR->sleft;              /* whether to search for boundary   */
  GEN->sright = PAR->sright;
  GEN->max_ivs = PAR->max_ivs;           /* maximum number of subintervals   */

  /* initialize variables */
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;
  GEN->dleft = -UNUR_INFINITY;
  GEN->dright = UNUR_INFINITY;
  GEN->Umax = 1.;
  GEN->iv = NULL;
  GEN->n_ivs = -1;        /* -1 indicates that there are no intervals at all */
  /*   GEN->guide_size = 0;  */
  /*   GEN->guide = NULL; */
  GEN->area = DISTR.area; /* we use the value in the distribution object as first guess */
  GEN->logPDFconstant = 0.;   /* rescaling constant for logPDF                  */


  GEN->aCDF = NULL;
  GEN->aTRX = NULL;
  GEN->Tmin = 0.;
  GEN->Tmax = 0.;


  /* allocate maximal array of intervals */
  /* [ Maybe we could move this into _unur_ptx_interval() ] */
  GEN->iv = _unur_xmalloc(GEN->max_ivs * sizeof(struct unur_ptx_interval) );

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_ptx_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_ptx_create() */

/*---------------------------------------------------------------------------*/

int
_unur_ptx_check_par( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check parameters of given distribution and method                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* points for searching computational domain */
  GEN->bleft = _unur_max(GEN->bleft_par,DISTR.domain[0]);
  GEN->bright = _unur_min(GEN->bright_par,DISTR.domain[1]);

  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

  /* domain of distribution (used when x with PDF(x)=0 are found) */
  GEN->dleft =  DISTR.domain[0];
  GEN->dright =  DISTR.domain[1];

  /* center of distribution */
  DISTR.center = unur_distr_cont_get_center(gen->distr);
  if (DISTR.center < GEN->dleft || DISTR.center > GEN->dright) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,
		"center must be in given domain of distribution");
    DISTR.center = _unur_max(DISTR.center,GEN->dleft);
    DISTR.center = _unur_min(DISTR.center,GEN->dright);
  }

  /* check center of distribution */
  if (gen->variant == PTX_VARIANT_PDF) {
    if (PDF(DISTR.center)<=0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "PDF(center) <= 0.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  return UNUR_SUCCESS;
} /* end of _unur_ptx_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_ptx_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE  ((struct unur_ptx_gen*)clone->datap)

  struct unur_gen *clone;
  int i;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_PTX_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* we do not need the table with CDF values; FIXME  */
  CLONE->aCDF = NULL;
  CLONE->aTRX = NULL;

  /* copy coefficients for Newton polynomial */
  CLONE->iv =  _unur_xmalloc((GEN->n_ivs+1) * sizeof(struct unur_ptx_interval) );
  memcpy( CLONE->iv, GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_ptx_interval) );

  for(i=0; i<=GEN->n_ivs; i++) {
    CLONE->iv[i].ti = _unur_xmalloc( GEN->order * sizeof(double) );
    CLONE->iv[i].zi = _unur_xmalloc( GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].ti, GEN->iv[i].ti, GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].zi, GEN->iv[i].zi, GEN->order * sizeof(double) );
  }

  /* copy guide table */
  /*   CLONE->guide = _unur_xmalloc( GEN->guide_size * sizeof(int) ); */
  /*   memcpy( CLONE->guide, GEN->guide, GEN->guide_size * sizeof(int) ); */

  return clone;

#undef CLONE
} /* end of _unur_ptx_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_ptx_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  int i;

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_PTX ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_PTX_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free guide table */
  /*   if (GEN->guide) free (GEN->guide); */

  /* free array for subintervals of adaptive Gauss-Lobatto integration */
  _unur_lobatto_free(&(GEN->aCDF));
  _unur_lobatto_free(&(GEN->aTRX));



  /* free tables of coefficients of interpolating polynomials */
  if (GEN->iv) {
    for(i=0; i<=GEN->n_ivs; i++){
      free(GEN->iv[i].ti);
      free(GEN->iv[i].zi);
    }
    free (GEN->iv);
  }

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_ptx_free() */

/*****************************************************************************/

/* int */
/* _unur_ptx_make_guide_table (struct unur_gen *gen) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* make a guide table for indexed search                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int i,j, imax; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PTX_GEN,UNUR_ERR_COOKIE); */

/*   /\* allocate blocks for guide table (if necessary). */
/*      (we allocate blocks for maximal guide table.) *\/ */
/*   GEN->guide_size = (int) (GEN->n_ivs * PTX_GUIDE_FACTOR); */
/*   if (GEN->guide_size <= 0) GEN->guide_size = 1; */
/*   GEN->guide = _unur_xrealloc( GEN->guide, GEN->guide_size * sizeof(int) ); */

/*   /\* maximum index for array of data *\/ */
/*   imax = GEN->n_ivs; */

/*   /\* create guide table *\/ */
/*   i = 0; */
/*   GEN->guide[0] = 0; */
/*   for( j=1; j<GEN->guide_size ;j++ ) { */
/*     while(GEN->iv[i+1].cdfi/GEN->Umax < j/(double)GEN->guide_size && i < imax) */
/*       i++; */
/*     if (i >= imax) break; */
/*     GEN->guide[j]=i; */
/*   } */

/*   /\* check i *\/ */
/*   i = _unur_min(i,imax); */

/*   /\* if there has been an round off error, we have to complete the guide table *\/ */
/*   for( ; j<GEN->guide_size ;j++ ) */
/*     GEN->guide[j] = i; */

/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */

/* } /\* end of _unur_ptx_make_guide_table() *\/ */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_eval_PDF (double x, struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* call to PDF.                                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x   ... argument of PDF                                            */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   PDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *distr = gen->distr;

  if (DISTR.logpdf != NULL) {
    return (exp((DISTR.logpdf)(x,distr) - GEN->logPDFconstant));
  }
  else {
    return ((DISTR.pdf)(x,distr));
  }

} /* end of _unur_ptx_eval_PDF() */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_eval_dTRX (double x, struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* call to FIXME.                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x   ... argument of PDF                                            */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   FIXME at x                                                           */
     /*----------------------------------------------------------------------*/
{
  double TRXx;
  double res,fx;

  /* derivative of transformation F_i^{-1}(F_o(x)): */
  /* (d/dx) (F_i^{-1}(F_o(x))) = f_o(x) / f_i(F_i^{-1}(F_o(x))) */

  TRXx = INVCDFIN(CDF(x));
  fx = _unur_cont_PDF(x,gen->distr);

  res = (fx / PDFIN(TRXx));

  printf("x=%g, ",x);
  printf("u=%g, ",CDF(x));
  printf("t=%g, ",TRXx);
  printf("fo(x)=%g, ",fx);
  printf("fi(t)=%g, ",PDFIN(TRXx));
  printf("res=%g\n",res);


  return res;

} /* end of _unur_ptx_eval_dTRX() */

/*---------------------------------------------------------------------------*/

