/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unuran.h                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and declares structures and function prototypes    *
 *         for all UNURAN methods                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         required for every application of UNURAN.                         *
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
#ifndef __UNURAN_H_SEEN
#define __UNURAN_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Basic header files                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* basic header files */
#include <stdio.h>
#include <stdlib.h>

/*---------------------------------------------------------------------------*/
/* compiler switches and defaults */
#include <unuran_config.h>

/*****************************************************************************/
/**  Typedefs                                                               **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* objects                                                                   */

typedef struct unur_distr UNUR_DISTR;    /* distribution object              */
typedef struct unur_par   UNUR_PAR;      /* parameters for generator         */
typedef struct unur_gen   UNUR_GEN;      /* generator object                 */

/*---------------------------------------------------------------------------*/
/* functions for continuous univariate PDF, CDF, and their derivatives       */

typedef double UNUR_FUNCT_CONT(double x, struct unur_distr *distr);
typedef double UNUR_FUNCT_DISCR(int x, struct unur_distr *distr);

/*---------------------------------------------------------------------------*/
/* functions for continuous multivariate PDF, CDF, and their gradients       */

typedef double UNUR_FUNCT_CVEC(double *x, struct unur_distr *distr);
typedef int UNUR_VFUNCT_CVEC(double *result, double *x, struct unur_distr *distr);

/*****************************************************************************/
/**  Declarations for uniform random number generators                      **/
/*****************************************************************************/

#include <x_urng.h>

/*****************************************************************************/
/**  Function prototypes for manipulating distribution objects              **/
/*****************************************************************************/

#include <distr.h>
#include <distr_cemp.h>
#include <distr_cont.h>
#include <distr_corder.h>
#include <distr_cvec.h>
#include <distr_cvemp.h>
#include <distr_discr.h>

/*****************************************************************************/
/**  Function prototypes for manipulating generator objects                 **/
/*****************************************************************************/

/* methods for discrete distributions */
#include <dari.h>
#include <dau.h>
#include <dgt.h>

/* methods for continuous distributions */
#include <arou.h>
#include <ninv.h>
#include <srou.h>
#include <ssr.h>
#include <tabl.h>
#include <tdr.h>
#include <utdr.h>
#include <empk.h>

/* methods for continuous multivariate distributions */
#include <vmt.h>
#include <vempk.h>

/* wrapper for special generators for standard distributions */
#include <cstd.h>     /* continuous */
#include <dstd.h>     /* discrete   */

/* wrapper for uniform random number generator */
#include <unif.h>


/*****************************************************************************/
/**  Invoke generators                                                      **/  
/*****************************************************************************/

#include <x_gen.h>

/*****************************************************************************/
/**  Distributions                                                          **/
/*****************************************************************************/

#include <unuran_distributions.h>

/*****************************************************************************/
/**  Debugging and Error messages                                           **/
/*****************************************************************************/

#include <x_errno.h>
#include <x_debug.h>

/*****************************************************************************/
/**  Additional header files for further function prototypes                **/
/*****************************************************************************/

#include <x_math.h>

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_H_SEEN */
/*---------------------------------------------------------------------------*/



