/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tdr_source.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros for method TDR                                    *
 *         (Transformed Density Rejection)                                   *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in tdr*.c files                                     *
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
#ifndef __TDR_SOURCE_H_SEEN
#define __TDR_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define TDR_VARMASK_T          0x000fu   /* indicates transformation         */
#define TDR_VAR_T_SQRT         0x0001u   /* T(x) = -1/sqrt(x)                */
#define TDR_VAR_T_LOG          0x0002u   /* T(x) = log(x)                    */
#define TDR_VAR_T_POW          0x0003u   /* T(x) = -x^c                      */

#define TDR_VARMASK_VARIANT    0x00f0u   /* indicates which variant          */
#define TDR_VARIANT_GW         0x0010u   /* original variant (Gilks&Wild)    */
#define TDR_VARIANT_PS         0x0020u   /* use proportional squeeze         */
#define TDR_VARIANT_IA         0x0030u   /* use immediate acceptance
					    (requires prop. squeeze)         */

#define TDR_VARFLAG_VERIFY     0x0100u   /* flag for verifying mode          */
#define TDR_VARFLAG_USECENTER  0x0200u   /* whether center is used as cpoint or not */
#define TDR_VARFLAG_USEMODE    0x0400u   /* whether mode is used as cpoint or not */
#define TDR_VARFLAG_PEDANTIC   0x0800u   /* whether pedantic checking is used */
#define TDR_VARFLAG_USEDARS    0x1000u   /* whether DARS is used in setup or not */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define TDR_DEBUG_IV           0x00000010u
#define TDR_DEBUG_SPLIT        0x00010000u
#define TDR_DEBUG_DARS         0x00020000u
#define TDR_DEBUG_SAMPLE       0x01000000u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define TDR_SET_CENTER         0x002u
#define TDR_SET_STP            0x004u
#define TDR_SET_N_STP          0x008u
#define TDR_SET_GUIDEFACTOR    0x010u
#define TDR_SET_C              0x020u
#define TDR_SET_MAX_SQHRATIO   0x040u
#define TDR_SET_MAX_IVS        0x080u
#define TDR_SET_USE_DARS       0x100u
#define TDR_SET_DARS_FACTOR    0x200u

/*---------------------------------------------------------------------------*/

#define GENTYPE "TDR"          /* type of generator                          */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.tdr         /* data for parameter object         */
#define GEN       gen->data.tdr         /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),&(gen->distr))   /* call to PDF         */
#define dPDF(x)   _unur_cont_dPDF((x),&(gen->distr))  /* call to derivative of PDF */

/*---------------------------------------------------------------------------*/
#endif   /* __TDR_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/








