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

/*****************************************************************************/
/**  Declarations for uniform random bumber generators                      **/
/*****************************************************************************/

#include <unur_urng.h>

/*****************************************************************************/
/**  Function prototypes for manipulating distribution objects              **/
/*****************************************************************************/

#include <unur_distr.h>

/*****************************************************************************/
/**  Function prototypes for manipulating generator objects                 **/
/*****************************************************************************/

/* methods for discrete distributions */
#include <unur_dau.h>
#include <unur_dis.h>

/* methods for continuous distributions */
#include <unur_arou.h>
#include <unur_ninv.h>
#include <unur_srou.h>
#include <unur_stdr.h>
#include <unur_tabl.h>
#include <unur_tdr.h>
#include <unur_utdr.h>

/* methods for continuous multivariate distributions */
#include <unur_rect.h>

/* wrapper for special generators for standard distributions */
#include <unur_cstd.h>     /* continuous */
#include <unur_dstd.h>     /* discrete   */

/* wrapper for uniform random number generator */
#include <unur_unif.h>


/*****************************************************************************/
/**  Invoke generators                                                      **/  
/*****************************************************************************/

UNUR_GEN *unur_init( UNUR_PAR *parameters );

int    unur_sample_discr(UNUR_GEN *generator);
double unur_sample_cont(UNUR_GEN *generator);
void   unur_sample_vec(UNUR_GEN *generator, double *vector);

void unur_free( UNUR_GEN *gen );

/*****************************************************************************/
/**  Additional header files for further function prototypes                **/
/*****************************************************************************/


#include <unur_misc.h>

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_H_SEEN */
/*---------------------------------------------------------------------------*/
