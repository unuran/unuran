/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_distr.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for handling               *
 *         distribution objects.                                             *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_unuran.h                                  *
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
/* indicate changed parameters                                               */

/* essential parameters */
#define UNUR_DISTR_SET_MASK_ESSENTIAL 0xffff0000

#define UNUR_DISTR_SET_DOMAIN         0x00010000
#define UNUR_DISTR_SET_STDDOMAIN      0x00020000 /* domain not truncated (for standard distributions) */

/* derived parameters */
#define UNUR_DISTR_SET_MASK_DERIVED   0x0000ffff

#define UNUR_DISTR_SET_MODE           0x00000001
#define UNUR_DISTR_SET_PDFAREA        0x00000002

/*---------------------------------------------------------------------------*/
/* call pdf's and cdf's                                                      */
/* (no checking for NULL pointer !)                                          */

#define _unur_cont_PDF(x,distr)   ((*((distr)->data.cont.pdf)) ((x),(distr)))
#define _unur_cont_dPDF(x,distr)  ((*((distr)->data.cont.dpdf))((x),(distr)))
#define _unur_cont_CDF(x,distr)   ((*((distr)->data.cont.cdf)) ((x),(distr)))

#define _unur_discr_PMF(x,distr)  ((*((distr)->data.discr.pmf))((x),(distr)))
#define _unur_discr_CDF(x,distr)  ((*((distr)->data.discr.cdf))((x),(distr)))

/*---------------------------------------------------------------------------*/
/* debuging routines for distributions                                       */

void _unur_distr_cont_debug( UNUR_DISTR *distribution, char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_corder_debug( UNUR_DISTR *order_statistics, char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_cemp_debug( UNUR_DISTR *distribution, char *genid, int printvector );
/* write info about distribution into logfile                                */

void _unur_distr_discr_debug( UNUR_DISTR *distribution, char *genid );
/* write info about distribution into logfile                                */

void _unur_distr_demp_debug( UNUR_DISTR *distribution, char *genid, int printvector );
/* write info about distribution into logfile                                */

/*---------------------------------------------------------------------------*/
/* check if parameter object is of correct type, return 0 otherwise       */

#define _unur_check_distr_object( distr,distrtype, rcode ) \
  do { \
    if ((distr)->type != UNUR_DISTR_##distrtype) { \
      _unur_warning((distr)->name,UNUR_ERR_DISTR_INVALID,""); \
      return rcode; } \
    COOKIE_CHECK(distr,CK_DISTR_##distrtype,rcode); } while (0)

/*---------------------------------------------------------------------------*/
