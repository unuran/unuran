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

  UNUR_FUNCT_CONT *pdf;         /* pointer to p.d.f.                         */
  UNUR_FUNCT_CONT *dpdf;        /* pointer to derivative of p.d.f.           */
  UNUR_FUNCT_CONT *cdf;         /* pointer to c.d.f.                         */

  double params[UNUR_DISTR_MAXPARAMS];  /* parameters of the p.d.f.          */
  int    n_params;              /* number of parameters of the pdf           */

  double norm_constant;         /* (log of) normalization constant for p.d.f.*/

  double mode;                  /* location of mode                          */
  double area;                  /* area below p.d.f.                         */
  double domain[2];             /* boundary of domain                        */

  int (*upd_mode)(struct unur_distr *distr);   /* funct for computing mode   */
  int (*upd_area)(struct unur_distr *distr);   /* funct for computing area   */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define object for multivariate continuous distribution                    */
struct unur_distr_cvec {

  UNUR_FUNCT_CVEC *pdf;         /* pointer to p.d.f.                         */
  UNUR_VFUNCT_CVEC *dpdf;       /* pointer to gradiant of p.d.f.             */

  double *params[UNUR_DISTR_MAXPARAMS]; /* parameters of the p.d.f.          */
  int    n_params;              /* number of parameters of the pdf           */

  double norm_constant;         /* (log of) normalization constant for p.d.f.*/

  double *mode;                 /* location of mode                          */
  double volume;                /* volume below p.d.f.                       */

  int (*upd_mode)(struct unur_distr *distr);   /* funct for computing mode   */
  int (*upd_volume)(struct unur_distr *distr); /* funct for computing volume */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define object for univariate discrete distribution                        */
struct unur_distr_discr {
  UNUR_FUNCT_DISCR *pmf;        /* pointer to probability mass function      */
  UNUR_FUNCT_DISCR *cdf;        /* pointer to c.d.f.                         */

  double params[UNUR_DISTR_MAXPARAMS];  /* parameters of the p.m.f.          */
  /* params[UNUR_DISTR_MAXPARAMS] is used to store normalization constants!! */
  int    n_params;              /* number of parameters of the pdf           */

  double norm_constant;         /* (log of) normalization constant for p.m.f.*/

  int domain[2];                /* boundary of domain                        */
  double area;                  /* area below p.m.f.                         */
  int (*upd_area)(struct unur_distr *distr);   /* funct for computing area   */

  int  (*init)(struct unur_par *par,struct unur_gen *gen);
                                /* pointer to special init routine           */
};

/*---------------------------------------------------------------------------*/
/* define object for empirical univariate discrete distribution              */
/* (given by probability vector)                                             */
struct unur_distr_demp {
  double *prob;                 /* pointer to probability vector             */
  int     n_prob;               /* length of probability vector              */
};

/*---------------------------------------------------------------------------*/
/* define object for empirical univariate constinuous distribution           */
/* (given by empirical sample)                                               */
struct unur_distr_cemp {
  double *sample;               /* pointer to sample                         */
  int     n_sample;             /* length of sample probability vector       */
};

/*---------------------------------------------------------------------------*/
/* define distribution object                                                */
struct unur_distr {
  union {             
    struct unur_distr_cont  cont;   /* univariate continuous distribution    */
    struct unur_distr_cvec  cvec;   /* multivariate continuous distribution  */
    struct unur_distr_discr discr;  /* univariate discrete distribution      */
    struct unur_distr_demp  demp;   /* empirical univariate discr. distr.    */
    struct unur_distr_cemp  cemp;   /* empirical univ. cont. distr. (sample) */
  } data;                           /* data for distribution                 */

  unsigned type;                    /* type of distribution                  */
  unsigned id;                      /* identifier for distribution           */
  const char *name;                 /* name of distribution                  */
  int dim;                          /* number of components of random vector */

  unsigned set;                     /* indicate changed parameters           */

  struct unur_distr *base;          /* pointer to distribution object for
				       derived distribution 
				       (e.g. order statistics)               */

#ifdef UNUR_COOKIES
  unsigned cookie;                  /* magic cookie                          */
#endif
};

/*---------------------------------------------------------------------------*/
