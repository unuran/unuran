/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen_source.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for handling               *
 *         generator objects.                                                *
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
/* Invoke generators (macros to avoid function calls)                        */  

#define _unur_init(par)               (par)->init(par)

#define _unur_sample_discr(gen)       (gen)->sample.discr(gen)
#define _unur_sample_cont(gen)        (gen)->sample.cont(gen)
#define _unur_sample_vec(gen,vector)  (gen)->sample.cvec(gen,vector)

#define _unur_free(gen)               do {if(gen) (gen)->destroy(gen);} while(0)

/*---------------------------------------------------------------------------*/
/* get type of transformation method                                         */

#define _unur_gen_is_discr(gen) ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ? 1 : 0 )
#define _unur_gen_is_cont(gen)  ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_CONT)  ? 1 : 0 )
#define _unur_gen_is_vec(gen)   ( (((gen)->method & UNUR_MASK_TYPE) == UNUR_METH_VEC)   ? 1 : 0 )

/*---------------------------------------------------------------------------*/
/* aux routine when no sampling routine is available                         */

double _unur_sample_cont_error( UNUR_GEN *gen );

/*---------------------------------------------------------------------------*/
/* copy (clone) generator objects                                            */

UNUR_GEN *_unur_arou_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_cstd_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_dari_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_dau_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_dgt_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_dsrou_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_dss_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_dstd_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_empk_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_empl_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_hinv_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_hrb_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_hrd_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_hri_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_mcorr_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_ninv_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_nrou_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_srou_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_ssr_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_tabl_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_tdr_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_unif_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_utdr_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_vempk_clone( const UNUR_GEN *gen );
UNUR_GEN *_unur_vmt_clone( const UNUR_GEN *gen );

/* no such routines:                                                         */
/* UNUR_GEN *_unur_auto_clone( const UNUR_GEN *gen );                        */

#define _unur_gen_clone(gen)    ((gen)->clone(gen))

/*---------------------------------------------------------------------------*/
