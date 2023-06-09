/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ptx_info.c                                                  *
 *                                                                           *
 *   Routines for creating info strings.                                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_ptx_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);

  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = %s\n",
		      (gen->variant==PTX_VARIANT_PDF) ? "PDF" : "CDF");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.trunc[0],DISTR.trunc[1]);
  _unur_string_append(info,"   center    = %g", unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]\n");
    else 
      _unur_string_append(info,"  [default]\n");
  }
  else {
    _unur_string_append(info,"\n");
  }

  if (help) {
    if ( !(distr->set & (UNUR_DISTR_SET_CENTER | UNUR_DISTR_SET_MODE )) ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
                          "You may provide a point near the mode as \"center\"."); 
  }
  _unur_string_append(info,"\n");
  

  /* method */
  _unur_string_append(info,"method: PTX (Polynomial interpolation based INVerse CDF)\n");
  _unur_string_append(info,"   order of polynomial = %d\n", GEN->order);
  switch (gen->variant) {
  case PTX_VARIANT_PDF:
    _unur_string_append(info,"   use PDF + Lobatto integration  %s\n",
			(gen->set & PTX_SET_VARIANT) ? "" : "[default]");
    break;
  case PTX_VARIANT_CDF:
    _unur_string_append(info,"   use CDF  %s\n",
			(gen->set & PTX_SET_VARIANT) ? "" : "[default]");
    break;
  }
  _unur_string_append(info,"\n");


  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   truncated domain = (%g,%g)\n",GEN->bleft,GEN->bright);
  /*   _unur_string_append(info,"   Prob(X<domain)   = %g\n", _unur_max(0,GEN->tailcutoff_left)); */
  /*   _unur_string_append(info,"   Prob(X>domain)   = %g\n", _unur_max(0,1.-GEN->tailcutoff_right)); */
  if (DISTR.cdf) {
    double max_error=1.; double MAE=1.;
    unur_ptx_estimate_error( gen, 10000, &max_error, &MAE );
    _unur_string_append(info,"   u-error         <= %g  (mean = %g)\n", max_error, MAE);
  }
  else {
    _unur_string_append(info,"   u-error            NA  [requires CDF]\n");
  }

  _unur_string_append(info,"   # intervals      = %d\n", GEN->n_ivs);
  _unur_string_append(info,"\n");
  

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   order = %d  %s\n", GEN->order,
 			(gen->set & PTX_SET_ORDER) ? "" : "[default]");

    _unur_string_append(info,"   u_resolution = %g  %s\n", GEN->u_resolution,
 			(gen->set & PTX_SET_U_RESOLUTION) ? "" : "[default]");
    
    _unur_string_append(info,"   boundary = (%g,%g)  %s\n", GEN->bleft_par, GEN->bright_par,
			(gen->set & PTX_SET_BOUNDARY) ? "" : "[default]");

    _unur_string_append(info,"   search for boundary: left=%s,  right=%s  %s\n",
			(GEN->sleft ? "TRUE":"FALSE"), (GEN->sright ? "TRUE":"FALSE"), 
			(gen->set & PTX_SET_BOUNDARY) ? "" : "[default]");

    _unur_string_append(info,"   maximum number of interval = %d  %s\n", GEN->max_ivs,
			(gen->set & PTX_SET_MAX_IVS) ? "" : "[default]");

    _unur_string_append(info,"\n");
  }


  /* Hints */
  if (help) {
    if ( GEN->order < MAX_ORDER )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase \"order\" to decrease #intervals");
    if (! (gen->set & PTX_SET_U_RESOLUTION) )
      _unur_string_append(info,"[ Hint: %s\n\t%s ]\n",
			  "You can decrease the u-error by decreasing \"u_resolution\".",
			  "(it is bounded by the machine epsilon, however.)");
    _unur_string_append(info,"\n");
  }

} /* end of _unur_tdr_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
