
/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

/*****************************************************************************/
/**  Sampling                                                               **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_sample_cvec( struct unur_gen *gen, double *rpoint )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  CONE *c;       /* cone for generating point */
  double gx;     /* distance of random point */
  double U;      /* uniformly distributed random number */
  double f, h;   /* value of density and hat at random point */
  int i,j;

  double *S = GEN->S;            /* working array for storing point on simples */

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  /* loop until random point is accepted */
  while( 1 ) { 

    /*.................................................................*/
    /** find a cone **/

    U = _unur_call_urng(gen->urng);    /* sample from uniform distribution */

    /* look up in guide table and search for cone */
    for( c = (GEN->guide)[(int) (U * GEN->guide_size)]; 
	 c!=NULL && c->Hsum < GEN->Htot * U;
	 c=c->next );

/* #if DEBUG > 0 */
/*     if( c == NULL ) { */
/*       /\* something is wrong. This should not happen *\/ */
/* /\*       WARNING( "c==NULL" ); *\/ */
/*       c=GEN->last_c;      /\* we simple use the last cone *\/ */
/*     } */
/*     if( DEBUG & DB_RPOINT ) */
/*       fprintf(LOG,"H = %g (%g) --> cone = %d",GEN->Htot * u,GEN->Htot,c->index); */
/* #endif */

    /*.................................................................*/
    /** get random point and distance of hyper plane **/

    /* get x value for marginal distribution of hat --> hyperplane */
    gx = unur_sample_cont(GEN_GAMMA) / (c->beta);
      
    /* nonnegative uniform random numbers with sum u_i = 1 */
    _unur_mvtdr_simplex_sample(gen, S);
      
    /* move point into center */
    for( i=0; i<GEN->dim; i++ ) rpoint[i] = GEN->center[i];
      
    /* calculate random point on chosen hyper-plane */
    for( j=0; j<GEN->dim; j++ ) {
      double x = gx * S[j] / c->gv[j];
      for( i=0; i<GEN->dim; i++ )
	rpoint[i] += x * (c->v[j])->coord[i];
    }

    /*.................................................................*/

/* #if RECTANGLE == 1 */
/*     /\* is point inside domain (rectangle) ? *\/ */
/* #endif */

    /*.................................................................*/
    /** accept or reject **/

    f = PDF(rpoint);                        /* density */
    h = T_inv( c->alpha - c->beta * gx );   /* hat */

    if( _unur_call_urng(gen->urng) * h <= f )
      return;
  }

} /* end of _unur_mvtdr_sample_cvec() */

/*-----------------------------------------------------------------*/

int
_unur_mvtdr_simplex_sample( const struct unur_gen *gen, double *U )
     /*----------------------------------------------------------------------*/
     /* sample point uniformly on standard simplex                           */
     /* point in standard simplex 0 <= u[0] <= ... <= u[dim-2] <= 1          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   U   ... array for storing result                                   */
     /*----------------------------------------------------------------------*/
{
  int dim = GEN->dim;

  if (dim == 2) {
    U[0] = _unur_call_urng(gen->urng);
    U[1] = 1. - U[0];
    return UNUR_SUCCESS;
  }
  /*................................................................*/
  if (dim == 3) {
    U[0] = _unur_call_urng(gen->urng);
    U[1] = _unur_call_urng(gen->urng);
    if( U[0] > U[1] ) {
      U[2] = U[0]; U[0] = U[1]; U[1] = U[2];
    }
    U[2] = 1. - U[1];
    U[1] = U[1] - U[0];
    return UNUR_SUCCESS;
  }
  /*................................................................*/
  if (dim >3) {
    int i,j;
    double U_aux;

    /* generate dim-1 numbers */
    for( i=0; i<dim-1; i++ )
      U[i] = _unur_call_urng(gen->urng);

    /* sort numbers (insertion sort) */
    for( i=1; i<dim-1; i++ ) {
      U_aux = U[i];
      for( j=i; j>0 && U[j-1] > U_aux; j-- )
	U[j] = U[j-1];
      U[j] = U_aux;
    }
    
    /* transform to affine coordinates */
    U[dim-1] = 1.;
    for( i=dim-1; i>0; i-- )
      U[i] -= U[i-1];

    return UNUR_SUCCESS;
  }

  /* else --> make error message!! */

  return UNUR_SUCCESS;
} /* end of _unur_mvtdr_simplex_sample() */

/*---------------------------------------------------------------------------*/
