/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for distributions                             *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
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
/* define object for univariate continuous distribution                      */
struct unur_distr_cont {

  _UNUR_FUNCTION_CONT *pdf;     /* pointer to p.d.f.                         */
  _UNUR_FUNCTION_CONT *dpdf;    /* pointer to derivative of p.d.f.           */
  _UNUR_FUNCTION_CONT *cdf;     /* pointer to c.d.f.                         */

  double params[UNUR_DISTR_MAXPARAMS + 1];  /* parameters of the p.d.f.      */
  /* params[UNUR_DISTR_MAXPARAMS] is used to store normalization constants!! */
  int    n_params;              /* number of parameters of the pdf           */

  double mode;                  /* location of mode                          */
  double area;                  /* area below p.d.f.                         */
  double domain[2];             /* boundary of domain                        */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define object for univariate discrete distribution                        */
struct unur_distr_discr {
  _UNUR_FUNCTION_DISCR *pmf;    /* pointer to probability mass function      */
  _UNUR_FUNCTION_DISCR *cdf;    /* pointer to c.d.f.                         */

  double params[UNUR_DISTR_MAXPARAMS + 1];  /* parameters of the p.d.f.      */
  /* params[UNUR_DISTR_MAXPARAMS] is used to store normalization constants!! */
  int    n_params;              /* number of parameters of the pdf           */

  double *prob;                 /* pointer to probability vector             */
  int     n_prob;               /* length of probability vector              */

  double domain[2];             /* boundary of domain                        */
  double area;                  /* area below p.d.f.                         */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define distribution object                                                */

struct unur_distr {
  union {             
    struct unur_distr_cont  cont;   /* univariate continuous distribution    */
    struct unur_distr_discr discr;  /* univariate discrete distribution      */
  } data;                           /* data for distribution                 */

  unsigned type;                    /* type of distribution                  */
  unsigned id;                      /* identifier for distribution           */
  const char *name;                 /* name of distribution                  */

  unsigned set;                     /* indicate changed parameters           */

#ifdef UNUR_COOKIES
  unsigned cookie;                  /* magic cookie                          */
#endif
};

/*---------------------------------------------------------------------------*/
