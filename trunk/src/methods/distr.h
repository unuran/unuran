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

/* 
   =DISTR  Distribution_objects  Handling distribution objects

   =UP TOP [30]

   =DESCRIPTION
      Objects of type @code{UNUR_DISTR} are used for handling
      distributions. All data about a distribution are stored on this
      object. UNURAN provides functions that return such objects for
      standard distributions 
      (@pxref{Stddist,Standard distributions,Standard distributions}).
      It is then possible to change this distribution object by
      various set calls. Moreover it is possible to build a
      distribution object entirely from scratch. For this purpose
      there exists an @command{unur_distr_<type>_new} call for each
      object type that returns an empty object of this type
      (eg. univariate contiuous) which can be filled with the
      appropriate set calls.
      
      Notice that there are essential data about a distribution,
      eg. the PDF, a list of (shape, scale, location) parameters for
      the distribution, and the domain of (the possibly truncated)
      distribution. 
      And there exist parameters that are/can be derived from these,
      eg. the mode of the distribution or the area below the given PDF 
      (which need not be normalized for many methods).
      UNURAN keeps track of parameters which are known. Thus if one of
      the essential parameters is changed all derived parameters are
      marked as unknown and must be set again if these are required
      for the chosen generation method.
      
      The library can handle truncated distributions, that is,
      distribution that are derived from (standard) distribution by
      simply restricting its domain to a subset. 
      However there is a subtle difference between changing the domain
      of a distribution object by a unur_distr_cont_set_domain() call
      and changing the (truncated) domain for an existing generator
      object. The domain of the distribution object is used to
      create the generator object with hats, squeezes, tables, etc.
      Whereas truncating the domain of an existing generator object
      need not necessarily require a recomputation of these data.
      Thus by a @command{unur_<method>_chg_truncated} call (if
      available) the sampling region is restricted to the subset of
      the domain of the given distribution object. However
      generation methods that require a recreation of the generator
      object when the domain is changed have a
      @command{unur_<method>_chg_domain} call instead. 
      For this call there are of course no restrictions on the
      given domain (i.e., it is possible to increase the domain of the
      distribution)
      (@pxref{Methods}, for details).

      For the objects provided by the UNURAN library of standard
      distributions, calls for updating these parameters exist (one for
      each parameter to avoid computational overhead since not all
      parameters are required for all generator methods).

      The calls listed below only handle distribution object. Since every
      generator object has its own copy of a distribution object, there are 
      calls for a chosen method that change this copy of distribution object.
      NEVER extract the distribution object out of the generator object and 
      run one of the below set calls on it.
      (How should the poor generator object know what has happend?)

   =EON
*/

/*

   =DISTR   AllDistr   Functions for all kinds of distribution objects

   =UP Distribution_objects [05]

   =END
*/

/*---------------------------------------------------------------------------*/
/* types of distribtuions                                                    */

enum {
  UNUR_DISTR_CONT  = 0x010u,     /* univariate continuous distribution       */ 
  UNUR_DISTR_CEMP  = 0x011u,     /* empirical univ. cont. distr. (a sample)  */ 
  UNUR_DISTR_CVEC  = 0x110u,     /* mulitvariate continuous distribution     */ 
  UNUR_DISTR_CVEMP = 0x111u,     /* empirical multiv. cont. distr. (sample)  */ 
  UNUR_DISTR_DISCR = 0x020u,     /* univariate discrete distribution         */ 
};

/*---------------------------------------------------------------------------*/

/*
  Parameters common to all distributions.
*/

/* =ROUTINES */

void unur_distr_free( UNUR_DISTR *distribution );
/* 
   Destroy a distribution object.
*/

int unur_distr_set_name( UNUR_DISTR *distribution, const char *name );
/* */

const char *unur_distr_get_name( UNUR_DISTR *distribution );
/* 
   Set and get name of distribution.
*/


int unur_distr_get_dim( UNUR_DISTR *distribution );
/* 
   Get number of components of random vector (its dimension).
   For univariate distributions the it returns @code{1} (of course).
*/


unsigned int unur_distr_get_type( UNUR_DISTR *distribution );
/* 
   Get type of @var{distribution}. 
   Possible types are
   @table @code
   @item UNUR_DISTR_CONT
   univariate continuous distributions
   @item UNUR_DISTR_CEMP
   empirical continuous univariate distributions (ie. samples)
   @item UNUR_DISTR_CVEC
   continuous mulitvariate distributions
   @item UNUR_DISTR_CVEMP
   empirical continuous multivaraiate distributions (ie. samples)
   @item UNUR_DISTR_DISCR
   discrete univariate distributions
   @end table

   Alternatively the @command{unur_distr_is_<TYPE>}
   calls can be used.
*/

int unur_distr_is_cont( UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is a continuous univariate distribution.
*/

int unur_distr_is_cvec( UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is a continuous multivariate distribution.
*/

int unur_distr_is_cemp( UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is an empirical continuous univariate distribution,
   i.e. a sample.
*/

int unur_distr_is_cvemp( UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is an empirical continuous multivariate
   distribution.
*/

int unur_distr_is_discr( UNUR_DISTR *distribution );
/* 
   TRUE if @var{distribution} is a discrete univariate distribution.
*/

/* =END */

/*---------------------------------------------------------------------------*/



