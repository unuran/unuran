/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl_sample.c                                                *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a unimodal distribution                                 *
 *      produce random variate X with its density                            *
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

double
_unur_tabl_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,INFINITY);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen->urng);

    /* look up in guide table and search for interval */
    iv =  GEN->guide[(int) (u * GEN->guide_size)];
    u *= GEN->Atotal;

    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,INFINITY);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->xmax <= iv->xmin) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u < iv->Asqueeze ) {
      /* below squeeze */
      return( iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze ); 
    }

    else {
      /* between spueeze and hat --> have to valuate PDF */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      fx = PDF(x);

      /* being above squeeze is bad. split interval. */
      if (GEN->n_ivs < GEN->max_ivs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  switch (_unur_tabl_split_interval( gen, iv, x, fx,(gen->variant & TABL_VARMASK_SPLIT)) ) {
	  case UNUR_SUCCESS:
	  case UNUR_ERR_SILENT:
	    /** TODO: it is not necessary to update the guide table every time. 
		But then (1) some additional bookkeeping is required and
		(2) the guide table method requires a acc./rej. step. **/
	    if (_unur_tabl_make_guide_table(gen) != UNUR_SUCCESS) {
	      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create guide table");
	      /* There is no chance to run out of this error! */
	      /* however, this should never happen.           */
	    }
	    break;
	  default:
	    /* condition for PDF is violated! */
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  }
	}
	else {
	  /* no more construction points (avoid to many second if statement above */
	  GEN->max_ivs = GEN->n_ivs;
	}
      }

      /* now accept or reject */
      u = _unur_call_urng(gen->urng);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	return x;
    }
  }

} /* end of _unur_tabl_sample() */

/*****************************************************************************/

double
_unur_tabl_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,INFINITY);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen->urng);

    /* look up in guide table and search for interval */
    iv =  GEN->guide[(int) (u * GEN->guide_size)];
    u *= GEN->Atotal;
    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,INFINITY);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->xmax <= iv->xmin) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u <= iv->Asqueeze ) {
      /* below squeeze */
      x = iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze;
      /* test whether PDF is monotone */
      fx = PDF(x);
      if (_unur_FP_greater(fx,iv->fmax))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. PDF not monotone in interval");
      if (_unur_FP_less(fx,iv->fmin))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. PDF not monotone in interval");
      /* at last return number */
      return x;
    }

    else {
      /* between spueeze and hat --> have to valuate PDF */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      fx = PDF(x);

      /* test whether PDF is monotone */
      if (_unur_FP_greater(fx,iv->fmax))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. PDF not monotone in interval");
      if (_unur_FP_less(fx,iv->fmin))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. PDF not monotone in interval");

      /* being above squeeze is bad. split interval. */
      if (GEN->n_ivs < GEN->max_ivs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  switch (_unur_tabl_split_interval( gen, iv, x, fx,(gen->variant & TABL_VARMASK_SPLIT)) ) {
	  case UNUR_SUCCESS:
	  case UNUR_ERR_SILENT:
	    /** TODO: it is not necessary to update the guide table every time. 
		But then (1) some additional bookkeeping is required and
		(2) the guide table method requires a acc./rej. step. **/
	    if (_unur_tabl_make_guide_table(gen) != UNUR_SUCCESS) {
	      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create guide table");
	      /* There is no chance to run out of this error! */
	      /* however, this should never happen.           */
	    }
	    break;
	  default:
	    /* condition for PDF is violated! */
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  }
	}
	else {
	  /* no more construction points (avoid to many second if statement above */
	  GEN->max_ivs = GEN->n_ivs;
	}
      }

      /* now accept or reject */
      u = _unur_call_urng(gen->urng);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	return x;
    }
  }

} /* end of _unur_tabl_sample_check() */

/*****************************************************************************/


