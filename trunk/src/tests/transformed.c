/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      transform.c                                                  *
 *                                                                           *
 *   sample from generator and transform to uniform by means of c.d.f.       *
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

#include <unur_tests.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

double
_unur_sample_cont_transformed( struct unur_gen *gen, 
			       double (*cdf)(double x, double *fparam, int n_fparam) )
/*---------------------------------------------------------------------------*/
/* sample from an univariate continuous generator and                        */
/* transform to uniform by means of the c.d.f.                               */
/*                                                                           */
/* parameters:                                                               */
/*   gen    ... pointer to generator object                                  */
/*   cdf    ... pointer to c.d.f. of distribution                            */
/*                                                                           */
/* return:                                                                   */
/*   uniform random number                                                   */
/*---------------------------------------------------------------------------*/
{
  double F, Fl, Fr, Fdelta; /* value of cdf at left and right boundary point */
  double *fparam;           /* list of parameters                            */
  int n_fparam;             /* number of parameters                          */

  /* read data from generator structure */
  switch(gen->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(gen,CK_AROU_GEN,0.);
    fparam = gen->data.arou.pdf_param;
    n_fparam = gen->data.arou.n_pdf_param;
    Fl = (gen->data.arou.bleft <= -INFINITY) ? 0. : cdf(gen->data.arou.bleft, fparam, n_fparam); 
    Fr = (gen->data.arou.bright >=  INFINITY) ? 1. : cdf(gen->data.arou.bright, fparam, n_fparam); 
    Fdelta = Fr - Fl;
    break;
  case UNUR_METH_CSTD:
    COOKIE_CHECK(gen,CK_CSTD_GEN,0.);
    fparam = gen->data.cstd.dist_param;
    n_fparam = gen->data.cstd.n_dist_param;
    Fl = 0.;
    Fr = 1.;
    Fdelta = 1.;
    break;
  case UNUR_METH_TABL:
    COOKIE_CHECK(gen,CK_TABL_GEN,0.);
    fparam = gen->data.tabl.pdf_param;
    n_fparam = gen->data.tabl.n_pdf_param;
    Fl = (gen->data.tabl.bleft <= -INFINITY) ? 0. : cdf(gen->data.tabl.bleft, fparam, n_fparam); 
    Fr = (gen->data.tabl.bright >=  INFINITY) ? 1. : cdf(gen->data.tabl.bright, fparam, n_fparam); 
    Fdelta = Fr - Fl;
    break;
  case UNUR_METH_TDR:
    COOKIE_CHECK(gen,CK_TDR_GEN,0.);
    fparam = gen->data.tdr.pdf_param;
    n_fparam = gen->data.tdr.n_pdf_param;
    Fl = (gen->data.tdr.bleft <= -INFINITY) ? 0. : cdf(gen->data.tdr.bleft, fparam, n_fparam); 
    Fr = (gen->data.tdr.bright >=  INFINITY) ? 1. : cdf(gen->data.tdr.bright, fparam, n_fparam); 
    Fdelta = Fr - Fl;
    break;
  case UNUR_METH_UTDR:
    COOKIE_CHECK(gen,CK_UTDR_GEN,0.);
    fparam = gen->data.utdr.pdf_param;
    n_fparam = gen->data.utdr.n_pdf_param;
    Fl = (gen->data.utdr.il <= -INFINITY) ? 0. : cdf(gen->data.utdr.il, fparam, n_fparam); 
    Fr = (gen->data.utdr.ir >=  INFINITY) ? 1. : cdf(gen->data.utdr.ir, fparam, n_fparam); 
    Fdelta = Fr - Fl;
    break;
  default:
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"unknown method. cannot compute c.d.f.");
    return -1.;
  }

  /* Fr - Fl <= 0. is a fatal error */
  if (Fdelta<=0.) {
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
    return -1.;
  }

  /* now run generator */
  F = cdf( unur_sample_cont(gen), fparam, n_fparam );

  /* return transformed uniform random number.
     rescale for trancated distributions */
  return ((F-Fl)/Fdelta);

} /* end of _unur_sample_cont_transformed() */

/*---------------------------------------------------------------------------*/






