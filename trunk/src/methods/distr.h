/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects         *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
/* Common variants for all special generators                                */

#define UNUR_STDGEN_INVERSION (~0u)

/*---------------------------------------------------------------------------*/
/* function prototypes for manipulating distribution object                  */

/*---------------------------------------------------------------------------*/
/* routines for all distribution objects                                     */

UNUR_DISTR *unur_distr_dup( UNUR_DISTR *distribution );
/* duplicate distribution object                                             */

void unur_distr_free( UNUR_DISTR *distribution );
/* destroy distribution object                                               */

int unur_distr_set_name( UNUR_DISTR *distribution, const char *name );
/* set name of distribution                                                  */

const char *unur_distr_get_name( UNUR_DISTR *distribution );
/* get name of distribution                                                  */

/*---------------------------------------------------------------------------*/
/* univariate continuous distributions                                       */

struct unur_distr *unur_distr_cont_new( void );
/* create a new distribution object                                          */

int unur_distr_cont_set_pdf( UNUR_DISTR *distribution, void *pdf );
/* set p.d.f. of distribution                                                */

double unur_distr_cont_pdf( UNUR_DISTR *distribution, double x );
/* evaluate p.d.f. of distribution at x                                      */

int unur_distr_cont_set_dpdf( UNUR_DISTR *distribution, void *dpdf );
/* set derivative of p.d.f. of distribution                                  */

int unur_distr_cont_set_cdf( UNUR_DISTR *distribution, void *cdf );
/* set c.d.f. of distribution                                                */

double unur_distr_cont_cdf( UNUR_DISTR *distribution, double x );
/* evaluate c.d.f. of distribution at x                                      */

int unur_distr_cont_set_params( UNUR_DISTR *distribution, double *params, int n_params );
/* set array of parameters for distribution                                  */

int unur_distr_cont_set_mode( UNUR_DISTR *distribution, double mode );
/* set mode of distribution                                                  */

double unur_distr_cont_get_mode( UNUR_DISTR *distribution );
/* get mode of distribution                                                  */

int unur_distr_cont_set_pdfarea( UNUR_DISTR *distribution, double area );
/* set area below p.d.f.                                                     */

int unur_distr_cont_set_domain( UNUR_DISTR *distribution, double left, double right );
/* set the left and right borders of the domain of the distribution          */

void _unur_distr_cont_debug( UNUR_DISTR *distribution, char *genid );
/* write info about distribution into logfile                                */

/*---------------------------------------------------------------------------*/
/* discrete univariate distributions                                         */

struct unur_distr *unur_distr_discr_new( void );
/* create a new distribution object                                          */

int unur_distr_discr_set_prob( UNUR_DISTR *distribution, double *prob, int n_prob );
/* set probability vector for distribution                                   */

void _unur_distr_discr_debug( UNUR_DISTR *distribution, char *genid, int printvector );
/* write info about distribution into logfile                                */

/*---------------------------------------------------------------------------*/
