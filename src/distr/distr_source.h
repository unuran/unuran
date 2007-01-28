/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_source.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for handling               *
 *         distribution objects.                                             *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_unuran.h                                  *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
#ifndef UNUR_DISTR_SOURCE_H_SEEN
#define UNUR_DISTR_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* indicate changed parameters                                               */

/* essential parameters */
#define UNUR_DISTR_SET_MASK_ESSENTIAL 0xffff0000u

#define UNUR_DISTR_SET_DOMAIN         0x00010000u
#define UNUR_DISTR_SET_DOMAINBOUNDED  0x00020000u /* domain is bounded */
#define UNUR_DISTR_SET_STDDOMAIN      0x00040000u /* domain not truncated (for standard distributions) */
#define UNUR_DISTR_SET_TRUNCATED      0x00080000u /* truncated distribution, i.e.
						    the domain of the distribution has been
						    restricted AFTER initializing
						    the generator object */

#define UNUR_DISTR_SET_MEAN           0x01000000u /* mean (vector for multivariate distr.) */
#define UNUR_DISTR_SET_COVAR          0x02000000u /* covariance matrix (for multivar. distr.) */
#define UNUR_DISTR_SET_COVAR_INV      0x04000000u /* inverse of covariance matrix (for multivar. distr.) */
#define UNUR_DISTR_SET_COVAR_IDENT    0x40000000u /* covariance matrix is identity matrix */
#define UNUR_DISTR_SET_CHOLESKY       0x08000000u /* cholesky factor of covariance matrix (for multivar. distr.) */
#define UNUR_DISTR_SET_RANKCORR       0x10000000u /* rank-correlation (for multivar. distr.) */
#define UNUR_DISTR_SET_RK_CHOLESKY    0x20000000u /* cholesky factor of covariance matrix (for multivar. distr.) */
#define UNUR_DISTR_SET_STDMARGINAL    0x00100000u /* standardized marginal distribution (for multivar. distr.) */
#define UNUR_DISTR_SET_MARGINAL       0x00200000u /* marginal distribution (for multivar. distr.) */

#define UNUR_DISTR_SET_GENERIC        0x00080000u /* generic parameter (can be used for any purpose) */


/* derived parameters */
#define UNUR_DISTR_SET_MASK_DERIVED   0x0000ffffu

#define UNUR_DISTR_SET_MODE           0x00000001u
#define UNUR_DISTR_SET_CENTER         0x00000002u
#define UNUR_DISTR_SET_PDFAREA        0x00000004u
#define UNUR_DISTR_SET_PMFSUM         0x00000008u
#define UNUR_DISTR_SET_PDFVOLUME      0x00000010u

/*---------------------------------------------------------------------------*/
/* call pdf's and cdf's                                                      */
/* (no checking for NULL pointer !)                                          */

#define _unur_cont_PDF(x,distr)     ((*((distr)->data.cont.pdf)) ((x),(distr)))
#define _unur_cont_dPDF(x,distr)    ((*((distr)->data.cont.dpdf))((x),(distr)))
#define _unur_cont_logPDF(x,distr)  ((*((distr)->data.cont.logpdf)) ((x),(distr)))
#define _unur_cont_dlogPDF(x,distr) ((*((distr)->data.cont.dlogpdf))((x),(distr)))
#define _unur_cont_CDF(x,distr)     ((*((distr)->data.cont.cdf)) ((x),(distr)))
#define _unur_cont_HR(x,distr)      ((*((distr)->data.cont.hr))  ((x),(distr)))

#define _unur_discr_PMF(x,distr)    ((*((distr)->data.discr.pmf))((x),(distr)))
#define _unur_discr_CDF(x,distr)    ((*((distr)->data.discr.cdf))((x),(distr)))

#define _unur_cvec_PDF(x,distr)        ((*((distr)->data.cvec.pdf)) ((x),(distr)))
#define _unur_cvec_dPDF(r,x,distr)     ((*((distr)->data.cvec.dpdf)) ((r),(x),(distr)))
#define _unur_cvec_pdPDF(x,c,distr)    ((*((distr)->data.cvec.pdpdf)) ((x),(c),(distr)))
#define _unur_cvec_logPDF(x,distr)     ((*((distr)->data.cvec.logpdf)) ((x),(distr)))
#define _unur_cvec_dlogPDF(r,x,distr)  ((*((distr)->data.cvec.dlogpdf)) ((r),(x),(distr)))
#define _unur_cvec_pdlogPDF(x,c,distr) ((*((distr)->data.cvec.pdlogpdf)) ((x),(c),(distr)))

/*---------------------------------------------------------------------------*/
/* check for existance of function pointers                                  */

#define _unur_cont_have_logPDF(distr)  (((distr)->data.cont.logpdf==NULL)?FALSE:TRUE)
#define _unur_cont_have_dlogPDF(distr) (((distr)->data.cont.dlogpdf==NULL)?FALSE:TRUE)

/*---------------------------------------------------------------------------*/
/* wrapper functions for PDF when only logPDF is given                       */

double _unur_distr_cont_eval_pdf_from_logpdf( double x, const struct unur_distr *distr );
double _unur_distr_cont_eval_dpdf_from_dlogpdf( double x, const struct unur_distr *distr );

double _unur_distr_cvec_eval_pdf_from_logpdf( const double *x, struct unur_distr *distr );
int _unur_distr_cvec_eval_dpdf_from_dlogpdf( double *result, const double *x, struct unur_distr *distr );
double _unur_distr_cvec_eval_pdpdf_from_pdlogpdf( const double *x, int coord, struct unur_distr *distr );

/*---------------------------------------------------------------------------*/
/* generic creator for distribution objects                                  */

struct unur_distr *_unur_distr_generic_new( void );

/*---------------------------------------------------------------------------*/
/* make clone of distribution objects                                        */

struct unur_distr *_unur_distr_cemp_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_cont_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_matr_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_cvec_clone ( const struct unur_distr *distr );
struct unur_distr *_unur_distr_cvemp_clone( const struct unur_distr *distr );
struct unur_distr *_unur_distr_discr_clone( const struct unur_distr *distr );

#define _unur_distr_clone(distr)    ((distr)->clone(distr))

/*---------------------------------------------------------------------------*/
/* destroy distribution object                                               */
#define _unur_distr_free(distr)    do {if (distr) (distr)->destroy(distr);} while(0)

/*---------------------------------------------------------------------------*/
/* debuging routines for distributions                                       */

void _unur_distr_cont_debug( const UNUR_DISTR *distribution, const char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_corder_debug( const UNUR_DISTR *order_statistics, const char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_cxtrans_debug( const UNUR_DISTR *cxtrans, const char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_cemp_debug( const UNUR_DISTR *distribution, const char *genid, unsigned printvector );
/* write info about distribution into logfile                                */

void _unur_distr_matr_debug( const UNUR_DISTR *distribution, const char *genid );
/* write info about matrix distribution into logfile                         */

void _unur_distr_cvec_debug( const UNUR_DISTR *distribution, const char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_condi_debug( const UNUR_DISTR *distribution, const char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_cvemp_debug( const UNUR_DISTR *distribution, const char *genid, unsigned printvector );
/* write info about distribution into logfile                                */

void _unur_distr_discr_debug( const UNUR_DISTR *distribution, const char *genid, unsigned printvector );
/* write info about distribution into logfile                                */

/*---------------------------------------------------------------------------*/
/* auxiliary routines                                                        */

/* test whether all marginals are equal or not  (returns TRUE or FALSE)      */
/* for dimesion 1, TRUE is returned.                                         */
/* WARNING: There is no checking of arguments in this function!              */
int _unur_distr_cvec_marginals_are_equal( struct unur_distr **marginals, int dim );

/*---------------------------------------------------------------------------*/
/* check if parameter object is of correct type, return 0 otherwise       */

#define _unur_check_distr_object( distr,distrtype, rcode ) \
  do { \
    if ((distr)->type != UNUR_DISTR_##distrtype) { \
      _unur_warning((distr)->name,UNUR_ERR_DISTR_INVALID,""); \
      return rcode; } \
    COOKIE_CHECK(distr,CK_DISTR_##distrtype,rcode); } while (0)

/*---------------------------------------------------------------------------*/
#endif   /* UNUR_DISTR_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
