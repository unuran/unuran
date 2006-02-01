
/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

/*****************************************************************************/
/**  Initialzation: Create Hat                                              **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mvtdr_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_MVTDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MVTDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mvtdr_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* we need an auxiliary generator for gamma random variates */
  GEN_GAMMA = _unur_mvtdr_gammagen( gen, (double)(GEN->dim) );
  if ( GEN_GAMMA == NULL ) {
      _unur_mvtdr_free(gen); return NULL;
  }

  /* make hat function */
  _unur_mvtdr_create_hat(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_mvtdr_debug_init(gen);
#endif

  /* free parameters */
  _unur_par_free(par);

  /* o.k. */
  return gen;

} /* end of _unur_mvtdr_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_mvtdr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MVTDR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_mvtdr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_MVTDR_GEN);

  /* dimension of distribution */
  GEN->dim = gen->distr->dim; 

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_mvtdr_sample_cvec;
  gen->destroy = _unur_mvtdr_free;
  gen->clone = _unur_mvtdr_clone;

  /* initialize counter and check given parameters */
  GEN->n_steps = 0;   /* no triangulation steps yet */
  GEN->steps_min = max( 0, PAR->steps_min );  /* minimum number of triangulation steps */
  /* check maximal number of cones */
  if ( (1 << (GEN->dim + GEN->steps_min)) > PAR->max_cones) {
    /*     WARNING( "number of cones raised to 2^(dim + T_STEPS_MIN)" ); */
    PAR->max_cones = 1 << (GEN->dim + GEN->steps_min);
  }

  /* initialize  pointers to lists */
  GEN->cone = NULL;
  GEN->last_cone = NULL;
  GEN->n_cone = 0;                      /* maximum number of vertices */

  GEN->vertex = NULL;
  GEN->last_vertex = NULL;
  GEN->n_vertex = 0;                      /* maximum number of vertices */

  GEN->etable = NULL;                   /* pointer to edge table */
  GEN->etable_size = 0;                    /* size of edge table */

  GEN->guide = NULL;
  GEN->guide_size = 0;
  
  /* initialize working arrays */
  GEN->S = _unur_xmalloc( GEN->dim * sizeof(double) );  /* for point on simples */
  GEN->g = _unur_xmalloc( GEN->dim * sizeof(double) );  /* vector g (direction of sweeping plane) */
  GEN->tp_coord = _unur_xmalloc( GEN->dim * sizeof(double) );  /* for coordinates of touching point of hat */
  GEN->tp_mcoord = _unur_xmalloc( GEN->dim * sizeof(double) );  /* for coordinates of touching point of hat moved into center */
  GEN->tp_Tgrad = _unur_xmalloc( GEN->dim * sizeof(double) );  /* for gradient of transformed density at tp */

  /* get center of the distribution */
  GEN->center = unur_distr_cvec_get_center(gen->distr);
 
  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_mvtdr_create() */

/*****************************************************************************/

struct unur_gen *
_unur_mvtdr_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_mvtdr_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MVTDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy data */
  CLONE->center = unur_distr_cvec_get_center(clone->distr);

  return clone;

#undef CLONE
} /* end of _unur_mvtdr_clone() */

/*****************************************************************************/

void
_unur_mvtdr_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  VERTEX *vt, *vt_next;
  CONE *c, *c_next;

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_MVTDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);


/*   db_close(gen->datap); */


  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* clear lists */

  /* hash table for edges */
  _unur_mvtdr_etable_free(gen);

  /* linked list of vertices */
  for (vt = GEN->vertex; vt != NULL; vt = vt_next) {
    vt_next = vt->next;
    free (vt->coord);      /* coordinates of vertex */
    free (vt);
  }

  /* linked list of cones */
  for (c = GEN->cone; c != NULL; c = c_next) {
      c_next = c->next;
      free (c->v);  /* list of vertices of the cone */
      free (c->center);   /* barycenter of cone */
      free (c->gv);  /* <g,v> for all vertices v */
      free (c);
  }

  /* guide table */
  if (GEN->guide) free (GEN->guide);

  /* working arrays */
  if (GEN->S) free (GEN->S); /* for point on simples */
  if (GEN->g) free (GEN->g);  /* vector g (direction of sweeping plane) */
  if (GEN->tp_coord) free (GEN->tp_coord);  /* for coordinates of touching point of hat */
  if (GEN->tp_mcoord) free (GEN->tp_mcoord);  /* for coordinates of touching point of hat moved into center */
  if (GEN->tp_Tgrad) free (GEN->tp_Tgrad);  /* for gradient of transformed density at tp */

  _unur_generic_free(gen);

} /* end of _unur_mvtdr_free() */

/*****************************************************************************/

struct unur_gen *
_unur_mvtdr_gammagen( struct unur_gen *gen, double alpha )
     /*----------------------------------------------------------------------*/
     /* create a gamma random variate generator with shape parameter alpha.  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to MVTDR generator object                        */
     /*   alpha ... shape parameter                                          */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *gammadistr;
  struct unur_par   *gammapar;
  struct unur_gen   *gammagen;
  double shape[1];

  shape[0] = alpha;
  gammadistr = unur_distr_gamma(shape,1);
  gammapar = unur_tdr_new( gammadistr );

  gammagen = unur_init( gammapar );
  _unur_distr_free( gammadistr );
  if (gammagen == NULL) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,
		"Cannot create aux Gamma generator");
    return NULL;
  }
  
  /* uniform random number generator and debugging flags */
  gammagen->urng = gen->urng;
  gammagen->debug = gen->debug;

  return gammagen;

} /* end of _unur_mvtdr_gammagen() */

/*****************************************************************************/
/*                                                                           */
/*   Hat.                                                                    */
/*                                                                           */
/*****************************************************************************/

int
_unur_mvtdr_create_hat( struct unur_gen *gen )
/* get hat function */
{
  int step;            /* triangulation steps */
  double Hi_bound;     /* lower bound on Hi for splitting cone */
  CONE *c;

  /* vertices of initial cones */
  _unur_mvtdr_initial_vertices(gen);

  /* initial cones */
  _unur_mvtdr_initial_cones(gen);

  /* execute minimal number of triangulation steps */
  for( step = 1; step <= GEN->steps_min; step++ )
    _unur_mvtdr_triangulate(gen,step,TRUE);

  /* compute optimal distance of touching points now */
  for( c = GEN->cone; c != NULL; c = c->next )
      _unur_mvtdr_tp_find (gen,c);

  /* some of cones with invalid hats (or too large volumes) must be split */
  while( _unur_mvtdr_triangulate(gen,step,FALSE) )
    step++;

  /* maximum number of triangulations yet */
  GEN->n_steps = step-1;

  /* compute cumulated volumes in all cones */
  GEN->Htot = 0.;                 /* accumulated sum of volumes */
  for( c=GEN->cone; c!=NULL; c=c->next ) {
    /* volume below hat */
    GEN->Htot += c->Hi;           /* volume below hat */
    c->Hsum = GEN->Htot;          /* accumulated sum of volumes */
  }

  /* bound for splitting cones */
  /* do until all cones have approx same hat volumes */
  Hi_bound = 1.5 * GEN->Htot / GEN->n_cone;

  /* and now check all the cones again */
  GEN->Htot = 0.;
  for( c=GEN->cone; c!=NULL; c=c->next ) {   /* all cones */
    while( Hi_bound < c->Hi ) { 
      /* we (must) split the cone again */
      _unur_mvtdr_cone_split(gen,c,c->level+1);
      /* and compute optimal touching point */
      _unur_mvtdr_tp_find (gen,c);
      _unur_mvtdr_tp_find (gen,GEN->last_cone);
    }
    GEN->Htot += c->Hi;           /* volume below hat */
    c->Hsum = GEN->Htot;          /* accumulated sum of volumes */
    if( c == GEN->last_cone ) break;
  }

  /* create guide table for finding cones */
  _unur_mvtdr_make_guide_table(gen);

  /* we do not need the hash table generated in triangulate() cone any more */
  if (GEN->dim > 2)
    _unur_mvtdr_etable_free(gen);

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_create_hat() */


/*****************************************************************************/
/*                                                                           */
/*   CONES.                                                                  */
/*                                                                           */
/*****************************************************************************/

int
_unur_mvtdr_initial_cones( struct unur_gen *gen )
/*
  get initial cones
*/
{
  int i,j,k;
  int max_c;       /* maximum number of cones */
  int dim = GEN->dim;
  CONE *c;
  VERTEX **ivtl;   /* list of initial vertices */
  VERTEX *vt;

  double logfak;
  /* factorial of (dim-1)*/
  logfak = 0.;
  for( i=2; i<dim; i++ ) 
    logfak += log((double) i);

  /* make array of initial vertices */
  ivtl = _unur_xmalloc( 2 * dim * sizeof(VERTEX*) );
  for (vt = GEN->vertex, i=0; i < 2*GEN->dim && vt!=NULL; vt = vt->next, i++)
    ivtl[i] = vt;    

  /* we have (at most) 2^dim initial cones */
  max_c = 1 << dim;

  /* we need vertices, index, volume for each cone */
  for( k=0; k<max_c; k++ ) {

    /* get new (empty) cone object */
    c = _unur_mvtdr_cone_new(gen);

    /* this is level 0 of triangulation */
    c->level = 0;

    /* each cone is incident to 'dim' edges.                           */
    /* The i-th edge of the cone is either GEN->v[i] or GEN->v[dim+i], */
    /* (the latter is equal to (-1)*(GEN->v[i])).                      */
    /* The indices of the vertices must be in ascending order.         */ 
    j = 0;
    for( i=0; i < dim; i++ )
      if (! ((k>>i)&1) * dim ) (c->v)[j++] = ivtl[i];
    for( i=0; i < dim && j < dim; i++ )
      if ( ((k>>i)&1) * dim ) (c->v)[j++] = ivtl[i + dim];

    /* determinant and volume of triangle [ = 1/((dim-1)!) ] */
    c->logdetf = -logfak;

    /* touching point not known yet */
    c->tp = -1.;   /* > 0. if and only if tp is computed !! */
  }

  /* free list of initial vertices */
  free (ivtl);

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_initial_cones() */

/*---------------------------------------------------------------------------*/

CONE *
_unur_mvtdr_cone_new( struct unur_gen *gen )
/*
  next allocated cone in list 
*/
{
  CONE *c; 

  /* allocate memory */
  if (GEN->cone == NULL) {
    c = GEN->last_cone = GEN->cone = _unur_xmalloc(sizeof(CONE));
  }
  else {
    c = GEN->last_cone = GEN->last_cone->next = _unur_xmalloc(sizeof(CONE));
  }
  c->next = NULL;

  /* list of vertices of the cone */
  c->v = _unur_xmalloc( GEN->dim * sizeof(VERTEX *) );

  /* barycenter of cone */
  c->center = _unur_xmalloc( GEN->dim * sizeof(double) );

  /* <g,v> for all vertices v */
  c->gv = _unur_xmalloc( GEN->dim * sizeof(double) );

  /* mark as invalid */
  c->tp = -1.;
  c->Hi = INFINITY;

  /* and update counter */
  ++(GEN->n_cone);

  /* return pointer to next vertex */
  return c;

} /* end of _unur_mvtdr_cone_new() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_cone_center( struct unur_gen *gen, CONE *c )
/*
  computer center of cone and normalize the corresponding vector it
*/
{
  int i,k;
  double norm;
  int dim = GEN->dim;
  
  norm = 0.;
  for( i=0; i<dim; i++ ) {
    c->center[i] = 0.;
    for( k=0; k<dim; k++ )
      c->center[i] += (c->v[k])->coord[i];        /* dim * barycenter */
    norm += c->center[i] * c->center[i];            /* norm ^2 */
  }
  
  /* norm --> 1 */
  norm = sqrt(norm);
  for( i=0; i<dim; i++ )
    c->center[i] /= norm;

  return UNUR_SUCCESS;
} /* end of _unur_mvtdr_cone_center() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_cone_params( struct unur_gen *gen, CONE *c )
/*
  compute parameters for hat for a cone
  (expect touching point and volume below hat)
*/
{
  double Tf,f;                  /* (transformed) density */
  double Tderf;                 /* T'(f(x)) */
  int i;                        /* aux variable */
  int dim = GEN->dim;           /* dimension */

  double *g = GEN->g;                  /* Vector g (direction of sweeping plane) */
  double *coord = GEN->tp_coord;              /* coordinates of touching point */
  double *mcoord = GEN->tp_mcoord;              /* coordinates of touching point moved into center */
  double *Tgrad = GEN->tp_Tgrad;              /* gradient of transformed density */

  /* coordinates of touching point */
  for( i=0; i<dim; i++ ) {
    coord[i] = c->tp * c->center[i];
    mcoord[i] = coord[i] + GEN->center[i];
  }

/* #if RECTANGLE == 1 */
/*     /\* check if point is in domain *\/ */
/*     for( i=0; i<N; i++ ) */
/*       if( coord[i] < GEN->rl[i] || coord[i] > GEN->ru[i] ) */
/* 	return TP_FZERO; */
/* #endif */

  /* density and its gradient */
  f = PDF(mcoord);
  dPDF(Tgrad,mcoord);

  /* check density */
  if( f < TOLERANCE )    /* f = 0. */
    return TP_FZERO;

  /* transformed density and its gradient */
  Tf = T(f);
  Tderf = T_deriv(f);
  for( i=0; i<dim; i++ )
    /* grad( T(f(x) ) = T'(f(x)) * grad(f(x)) */    
    Tgrad[i] *= Tderf;

  /* parameters alpha and beta */
  c->alpha = Tf - _unur_vector_scalar_product(dim,Tgrad,coord);
  c->beta = _unur_vector_norm(GEN->dim,Tgrad);

  /* |Tgrad| must not be to small */
  if( c->beta < TOLERANCE )
    return TP_GZERO;

  /* compute g and <g,v> for all vertices of cone */
  /* vector g = - grad(T(f)) / |grad(T(f))| */
  for( i=0; i<dim; i++ )
    g[i] = - Tgrad[i] / c->beta;

  /* <g,v> for each vertex v of cone and */
  /* parameter a1 for volume of cone */
  c->logai = c->logdetf;
  for( i=0; i<dim; i++ ) {
    (c->gv)[i] = _unur_vector_scalar_product(dim,g,(c->v[i])->coord);   /* <g,v> */
    if( c->gv[i] < TOLERANCE )
      return TP_UGRAD;
    else
      c->logai -= log(c->gv[i]);
  }

/* #if RECTANGLE == 1 */
/*   /\* at last calculate height of pyramid *\/ */
/* #if FIND_TP_FUNCTION == 1 */
/*   get_height(c);  /\* this is expensive for calculation for every touching point !!! *\/ */
/*   /\* TODO: approximate gat_height with   max_{vertices of rectangle} || vertex - mode || *\/ */
/* #endif */
/* #endif */

  /* return error code */
  return TP_VALID;

} /* end of _unur_mvtdr_cone_params() */

/*---------------------------------------------------------------------------*/

double _unur_mvtdr_cone_logH( struct unur_gen *gen, CONE *c )
     /* calculate log of volume below hat for given touching point */
     /* volume below hat for cone c constructed at touching point located at t.
       t > 0.
     */
{
  double logH;

  /* compute parameters for cone */
  if( _unur_mvtdr_cone_params(gen,c) != TP_VALID )
    /* something is wrong: beta = 0 and/or <g,v> <= 0 */
    return INFINITY;

  /* compute log of volume below hat */
  logH = c->alpha - GEN->dim * log(c->beta) + c->logai;

  /* check for numerical errors (alpha or beta too small) */
  if (_unur_isfinite(logH))
    return logH;

  else {
    /* we mark this case by setting tp = -1 */
    c->tp = -1.;
    return INFINITY;
  }
} /* end of _unur_mvtdr_cone_logH() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_cone_split( struct unur_gen *gen, CONE *c, int step )
/*
  split a cone along "oldest" edge
*/
{
  CONE *newc;       /* new cone */
  VERTEX *newv;      /* new vertex */
  int dim = GEN->dim;  /* dimension */
  int i;

  if (dim == 2)
    /* there is only one edge in this cone */
    newv = _unur_mvtdr_vertex_on_edge(gen,c->v);
  else
    /* find "oldest" edge, read center of edge from table (or computer) */
    newv = _unur_mvtdr_etable_find_or_insert(gen,c->v);

  /* construct two new cones */

  /* first cone */
  newc = _unur_mvtdr_cone_new(gen);                             /* pointer to new cone */
  newc->level = step;                             /* level */
  for (i=0; i<dim-1; i++)
    newc->v[i] = c->v[i+1]; /* copy list of vertices to new cone */
  newc->v[dim-1] = newv;                          /* the new vertex */
  newc->logdetf = c->logdetf - log(2.*newv->norm);        /* the log of determinant of spanning vectors */
  newc->tp = c->tp;                               /* distance of touching point remains unchanged */

  /* second cone */
  c->level = step;                                /* level */
  for (i=0; i<dim-2; i++)
    c->v[i+1] = c->v[i+2];             /* shift list of vertices */
  (c->v)[dim-1] = newv;                             /* the new vertex */
  c->logdetf = newc->logdetf;                           /* the determinant */

  /* store maximal triangulation level for debugging */
  GEN->n_steps = max(GEN->n_steps, step); 

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_cone_split() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_triangulate( struct unur_gen *gen, int step, int all )
/*
  make one triangulation step
*/
{
  int k,nc;
  CONE *c;
  int dim = GEN->dim;  /* dimension */

  if (dim > 2) {
    /* We need a hash table for storing the edges.                                 */
    /* size of table = maximum number of vertices in current triangulation cycle.  */
    /* only required at begining of new tiangulation iteration, i.e. step = dim-1. */
    /* length of cycle dim-1.                                                      */
    /* only necessary if dim > 2                                                   */
    if( step % (dim-1) == 1 )
      /* TODO: Warum wird das hier IMMER ausgefuehrt ? */
      _unur_mvtdr_etable_new(gen, _unur_mvtdr_number_vertices(gen, (step/(dim-1)+1)*(dim-1) ));
  }

  /*   number of cones before triangulation */
  nc = GEN->n_cone;

  /*   triangulate every cone */
  for( k=0, c=GEN->cone; k<nc; k++ ) {
    if( all )
      _unur_mvtdr_cone_split(gen,c,step);         /* split cone */
    else if ( c->tp < 0. ) {
      _unur_mvtdr_cone_split(gen,c,step);       /* split cone */
      _unur_mvtdr_tp_find (gen,c);
      _unur_mvtdr_tp_find (gen,GEN->last_cone);
    }
    /* next cone */
    c=c->next;
  }

  /* return number of new cones */
  return (GEN->n_cone - nc);

} /* end of _unur_mvtdr_triangulate() */


/*****************************************************************************/
/*                                                                           */
/*   Optimal distance for touching points                                    */
/*                                                                           */
/*****************************************************************************/

double _unur_mvtdr_tp_min (double t, void *p )
     /* 
	function to be minimized 
     */
{
  /* unpack arguments */
  TP_ARG *a = p;
  /* set new construction point for hat */
  (a->c)->tp = a->t = t;
  /* compute volume below hat in cone and return result */
  a->logH = _unur_mvtdr_cone_logH (a->gen, a->c);    /* max or min ?? */
  /* result o.k. ? */
  a->valid = (_unur_isfinite(a->logH)) ? TRUE : FALSE;

  if( !(a->valid) )
    /* we mark this case by setting tp = -1 */
    (a->c)->tp = -1.;
  
  return a->logH;
} /* end of _unur_mvtdr_tp_min() */

double _unur_mvtdr_tp_min_aux(double t, void *p) { 
  return (- _unur_mvtdr_tp_min(t, p) );
}

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_tp_find( struct unur_gen *gen, CONE *c )
/*
  find optimal touching point for cone
  Brent's method
*/
{
  struct unur_funct_generic tpaux;
  TP_ARG a[3];    /* left, middle and right point of bracket for Brent algorithm */
  int i;

  /* compute center */
  _unur_mvtdr_cone_center(gen,c);

  /* pack arguments for finding minimum */
  for (i=0; i<3; i++) { a[i].c = c; a[i].gen = gen; }

  /* find proper touching point */
  if( ! _unur_mvtdr_tp_search(gen,a) )
    /* no proper point found */
    return TP_EMPTY;

  /* find "bracket" for Brent's algorithm */
  switch( _unur_mvtdr_tp_bracket(gen,a) ) {      /* searching for intervall that contains minimum */
  case TP_BRACKET:                 /* bracket found */
    /* make auxiliary function for Brent's algorithms */
    tpaux.f = _unur_mvtdr_tp_min_aux;
    tpaux.params = a+1;
    c->tp = _unur_util_brent( tpaux, a[0].t, a[2].t, a[1].t, FIND_TP_TOL);
    c->Hi = exp(a[1].logH);
    /* TODO: error checking !! */
    break;                         /* c->tp already set by tp_min() */
  case TP_LEFT:                    /* minimum in left point */
    c->tp = a[0].t;
    c->Hi = exp(a[0].logH);
    break;
  case TP_MIDDLE:                  /* minimum in middle point */
    _unur_mvtdr_tp_min(a[1].t, a+1);  /* calculate volume function */
    c->Hi = exp(a[1].logH);
    break;
  case TP_RIGHT:                   /* minimum in right point */
    c->tp = a[2].t;
    c->Hi = exp(a[2].logH);
    break;
  case TP_EMPTY:                   /* no proper touching point found */
  default:
    c->tp = -1.;
    return UNUR_FAILURE;
  }

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_tp_find() */

/*-----------------------------------------------------------------*/

int
_unur_mvtdr_tp_search( struct unur_gen *gen, TP_ARG *a )
/* 
   search for proper touching point 
   returns 0 if no proper point could be found
   1 otherwise
   a[0], a[1], a[2] ... left, middle and right point of intervall
   the proper touching point is stored on a[1].
   (see find_bracket() )
*/
{
  int i;      /* aux counter */

/* #if RECTANGLE == 1 */
/*   /\* TODO: (a+1)->t has to be set to min of  FIND_TP_START and c->height if known *\/ */
/* #endif */

  /** search from 0 --> infinity **/

  /* initialize boundary of intervall */
  a[0].t = 0.;                  /* x[0] must >= 0. */
  a[1].t = FIND_TP_START;       /* starting point for searching proper touching point */
  a[2].t = -1.;                 /* not known. marked by setting to -1. */

  for( i=1; i <= max(1,FIND_TP_STEPS_R); i++ ) {
    /* TODO: if vol==0 for a point we need not continue */
    _unur_mvtdr_tp_min(a[1].t, a+1);  /* calculate volume function */

    if( ! a[1].valid ) {       /* check validity of touching point */
      /* not a proper touching point */
      a[0].t = a[1].t;
      a[1].t *= FIND_TP_STEP_SIZE;
    }
    else
      return 1;
  }
      
  /** search from infinity --> 0 **/

  /* initialize boundary of intervall */
  a[0].t = 0.;                  /* x[0] must >= 0. */
  a[1].t = FIND_TP_START / FIND_TP_STEP_SIZE;/* starting point for searching proper touching point */
  a[2].t = FIND_TP_START;       /* t[2] must >= t[1]. */

  for( i=0; i <= max(0,FIND_TP_STEPS_L); i++ ) {
    /* TODO: if vol==0 for a point we need not continue */
    _unur_mvtdr_tp_min(a[1].t, a+1);  /* calculate volume function */

    if( ! a[1].valid ) {       /* check validity of touching point */
      /* not a proper touching point */
      a[2].t = a[1].t;
      a[1].t /= FIND_TP_STEP_SIZE;
    }
    else
      return 1;
  }

  /* no proper touching point found */
  return 0;

} /* end of _unur_mvtdr_tp_search() */

/*-----------------------------------------------------------------*/

int 
_unur_mvtdr_tp_bracket( struct unur_gen *gen, TP_ARG *a )
/* 
   search for proper bracket of minimum of function tp_min()
   returns 0 if no proper bracket could be found
   1 otherwise
   a[0], a[1], a[2] ... left, middle and right point of intervall
*/
{
  int i;                 /* aux variable */
  double tleft, tright;  /* left boundary of searching region */

  /** left point of intervall **/

  /* initialize boundary of searching region and set starting point */
  tleft = a[0].t;
  a[0].t = a[1].t / 2.;

  /* search */
  for( i=1; i <= max(1,FIND_TP_STEPS_LB); i++ ) {
    _unur_mvtdr_tp_min(a[0].t, a);  /* calculate volume function */

    if( ! a[0].valid ) {
      /* a[0] not a proper touching point */
      tleft = a[0].t;                      /* change boundary of searching region */
      a[0].t += (a[1].t - a[0].t) / 2.;   /* try another one */
    }

    else if( a[0].logH <= a[1].logH ) {
      /* a[0] is proper touching point, but ... */
      a[2].t = a[1].t; a[2].logH = a[1].logH; a[2].valid = TRUE;
      a[1].t = a[0].t; a[1].logH = a[0].logH; a[1].valid = TRUE;
      a[0].t = tleft + (a[0].t - tleft)*0.5;
    }
    else  /* all right: a[0].logH > a[1].logH */
      break;
  }

  /* search successful ? */
  if( ! a[0].valid )
    /* no proper touching point on left side */
    return TP_MIDDLE;   /* ????? */
  if( a[0].logH <= a[1].logH )
    /* vol(left) <= vol(middle) */
    return TP_LEFT;
  
  /** right point of intervall **/

  /* initialize a[2] if necessary */
  if( a[2].t < 0. )
    a[2].t = 1.1 * a[1].t;
  tright = -1.;    /* no right boundary known yet */
  tleft = a[1].t;

  /* search */
  for( i=1; i <= max(1,FIND_TP_STEPS_UB); i++ ) {
    _unur_mvtdr_tp_min(a[2].t, a+2);  /* calculate volume function */
    if( ! a[2].valid ) {
      /* a[2] not a proper touching point */
      tright = a[2].t;
      a[2].t = (tleft + a[2].t) * 0.5;   /* try another one */
    }
    else if( a[2].logH <= a[1].logH ) {
      /* move right */
      tleft = a[2].t;
      a[2].t = (tright < 0.) ? a[2].t * 2. : (tright + a[2].t) * 0.5;
    }
    else    /* all right: vol(right) > vol(middle) */
      break;
  }

  /* search successful ? */
  if( ! a[2].valid )
    /* no proper touching point on right side */
    return TP_MIDDLE; 
  if( a[2].logH <= a[1].logH )
    /* f(right) <= f(middle) */
    return TP_RIGHT;

  /* we have found a bracket */
  return TP_BRACKET;

} /* end of _unur_mvtdr_tp_bracket() */

/*****************************************************************************/
/*                                                                           */
/*   VERTICES.                                                               */
/*                                                                           */
/*****************************************************************************/

int
_unur_mvtdr_initial_vertices( struct unur_gen *gen )
/*
  get vertices of initial cones and return number of vertices
*/
{
  VERTEX *vt;
  int i,k,d;

  /* unit vectors e_k  and -e_k */
  for ( d=1.; d > -2; d -= 2.) {
    for( k=0; k<GEN->dim; k++ ) {
      vt = _unur_mvtdr_vertex_new(gen);
      for( i=0; i<GEN->dim; i++ ) {
	/* coordinates */
	(vt->coord)[i] = (i==k) ? d : 0.;
      }
      /* all vectors have norm 1. */
      vt->norm = 1.;
    }
  }

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_initial_vertices() */

/*---------------------------------------------------------------------------*/

VERTEX *
_unur_mvtdr_vertex_new( struct unur_gen *gen )
/*
  next allocated vertex in list 
*/
{
  /* allocate memory */
  if (GEN->vertex == NULL) {
    GEN->last_vertex = GEN->vertex = _unur_xmalloc(sizeof(VERTEX));
  }
  else {
    GEN->last_vertex = GEN->last_vertex->next = _unur_xmalloc(sizeof(VERTEX));
  }
  (GEN->last_vertex)->next = NULL;
  
  /* coordinates of vertex */
  (GEN->last_vertex)->coord = _unur_xmalloc(GEN->dim * sizeof(double));

  /* index of vertex */
  (GEN->last_vertex)->index = GEN->n_vertex;
  /* and update counter */
  ++(GEN->n_vertex);

  /* return pointer to next vertex */
  return GEN->last_vertex;

} /* end of _unur_mvtdr_vertex_new() */

/*---------------------------------------------------------------------------*/

VERTEX *
_unur_mvtdr_vertex_on_edge( struct unur_gen *gen, VERTEX **vl )
/*
  compute new vertex on edge and return its index
*/
{
  int i;
  VERTEX *newv;          /* pointer to new vertex */

  /* get an empty vertex */
  newv = _unur_mvtdr_vertex_new(gen);

  /* barycenter of edge */
  for( i=0; i<GEN->dim; i++ )
    newv->coord[i] =
      0.5 * ( ((vl[0])->coord)[i] + ((vl[1])->coord)[i] );

  /* norm */
  newv->norm = _unur_vector_norm(GEN->dim, newv->coord);

  /* norm --> 1 */
  for( i=0; i<GEN->dim; i++ )
    newv->coord[i] /= newv->norm;

  return newv;

} /* end of _unur_mvtdr_vertex_on_edge() */

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_number_vertices( struct unur_gen *gen, int level )
/*
  calculate (approximate) number of vertices in given triangulation level
  (WARNING! The numbers are found by computer experiments)
*/
{
  if (level < 0 || GEN->dim < 2) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return -1;  /** TODO: errorcode **/
  }

  switch (GEN->dim) {

  case 2: {
    return (1 << (level+2));
  }
  case 3: {
    static int nv[]={ 6, 13, 18, 40, 66,142,258,538,1026,2098,4098,8290,16386,32962,65538,131458,262146};
    return nv[min(level,16)];
  }
  case 4: {
    static int nv[]={ 8, 19, 25, 32, 80,128,192,456, 824,1408,3120,5968,11008,23264,45600, 87552};
    return nv[min(level,15)];
  }
  case 5: {
    static int nv[]={10, 26, 33, 41, 50,140,220,321, 450,1186,2158,3636, 5890,13970,27130};
    return nv[min(level,14)];
  }
  case 6: {
    static int nv[]={12, 34, 42, 51, 61, 72,224,348, 501, 681, 912,2660, 4896, 8254};
    return nv[min(level,13)];
  }
  case 7: {
    static int nv[]={14, 43, 52, 62, 73, 85, 98,336, 518, 743, 985,1289, 1666};
    return nv[min(level,12)];
  }
  case 8: {
    static int nv[]={16, 53, 63, 74, 86, 99,113,128, 480, 736,1059};
    return nv[min(level,10)];
  }
  case 9: {
    static int nv[]={18, 64, 75, 87,100,114,129,145, 162, 660};
    return nv[min(level,9)];
  }
  case 10: {
    static int nv[]={20, 76, 88,101,115,130,146,163, 181, 200};
    return nv[min(level,9)];
  }
  case 11: {
    static int nv[]={22, 89,102,116,131,147,164,182, 201, 221, 242};
    return nv[min(level,10)];
  }
  default: { /* dim >= 12 */
    static int nv[]={24,103,117,132,148,165,183,202, 222, 243, 265, 288};
    return nv[min(level,11)];
  }
  }

} /* end of _unur_mvtdr_number_vertices() */


/*****************************************************************************/
/*                                                                           */
/*   hash table for storing EDGES.                                           */
/*                                                                           */
/*****************************************************************************/

/* hash function for edge table */
#define _unur_mvtdr_etable_hash(x,y)  ( (3*((x)+(y))/2) % GEN->etable_size )

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_etable_new( struct unur_gen *gen, int size )
/*
  make new hash table
*/
{
  int n;

  /* first clear hash table (if necessary) */
  _unur_mvtdr_etable_free(gen);

  /* set size of edge table */
  GEN->etable_size = size;

  /* make root */
  GEN->etable = _unur_xmalloc( size * sizeof(E_TABLE*) );

  /* initialize table */
  for (n = 0; n< size; n++) 
    GEN->etable[n] = NULL;
  
  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_etable_new() */

/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_etable_free( struct unur_gen *gen )
/*
  free hash table
*/
{
  int i;
  E_TABLE *et, *et_next;

  if( GEN->etable == NULL )
    /* nothing to do */
    return;

  /* clear all branches */
  for( i=0; i<GEN->etable_size; i++ ) {
    /* free branch in table (linked list) */
    for (et = GEN->etable[i]; et != NULL; et = et_next) {
      et_next = et->next;
      free (et);
    }
  }

  /* clear root */
  free( GEN->etable );
  GEN->etable = NULL;
  GEN->etable_size = 0;

} /* end if _unur_mvtdr_etable_free() */

/*---------------------------------------------------------------------------*/

VERTEX *
_unur_mvtdr_etable_find_or_insert( struct unur_gen *gen, VERTEX **vidx )
/*
  find or insert entry in hash table
  it returns pointer to corresponding vertex (barycenter of edge)
  we use "oldest" edge which is given by vertices vidx[0] and vidx[1]
*/
{
  E_TABLE *pet, *pet_last;  /* table entry */
  int idx[2];               /* store indices */
  int hidx;                 /* hash number */

  /* check pointer */
  CHECK_NULL(GEN->etable,NULL);

  /* hash number */
  idx[0] = vidx[0]->index;
  idx[1] = vidx[1]->index;
  hidx = _unur_mvtdr_etable_hash(idx[0],idx[1]);

  /* get branch of hash table */
  pet = pet_last = *(GEN->etable + hidx);

  /* now find entry in branch */
  while( pet != NULL ) {
    if( pet->index[0] == idx[0] && pet->index[1] == idx[1] )
      break;   /* found ! */
    pet_last = pet;
    pet =  pet->next;
  }

  if( pet == NULL ) {
    /* we have not found the index */

    /* new entry in hash table */
    pet = _unur_xmalloc( sizeof(E_TABLE) );
    pet->next = NULL;
    if (pet_last == NULL)
      *(GEN->etable + hidx) = pet;
    else
      pet_last->next = pet;

    /* insert data of new edge */
    /* indices of incident vertices */
    pet->index[0] = idx[0];
    pet->index[1] = idx[1];

    /* compute new vertex */
    pet->vertex = _unur_mvtdr_vertex_on_edge(gen,vidx);
  }

  /* return pointer to (new) edge */
  return pet->vertex;

} /* end of _unur_mvtdr_etable_find_or_insert() */

/*****************************************************************************/

int
_unur_mvtdr_make_guide_table( struct unur_gen *gen )
/* 
   create guide table 
*/
{
  int j;
  CONE *c;

  /* memory for the guide table */
  GEN->guide_size = GEN->n_cone * GUIDE_TABLE_SIZE;
  GEN->guide = _unur_xmalloc (GEN->guide_size * sizeof(CONE*));
  /* initialize table */
  for( j = 0; j < GEN->guide_size ; j++ )
    GEN->guide[j] = NULL;

  /* make table */
  for( c=GEN->cone, j=0; c!=NULL && j<GEN->guide_size; j++ ) {
    while( c->Hsum / GEN->Htot < (double) j / GEN->guide_size )
      c=c->next;
    (GEN->guide)[j] = c;
    if( c == GEN->last_cone ) break;
  }

  /* is there an error ? */
  if( j<GEN->guide_size )
    /* this should not happen */
    for( ; j<GEN->guide_size; j++ )
      (GEN->guide)[j] = GEN->last_cone;

  return UNUR_SUCCESS;

} /* end of _unur_mvtdr_make_guide_table() */

/*---------------------------------------------------------------------------*/
