/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unuran_errno.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines error codes.                                              *
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
#ifndef __UNURAN_ERRNO_H_SEEN
#define __UNURAN_ERRNO_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* error types                                                               */

enum { 

  /** distribution object **/
  UNUR_ERR_DISTR_SET      = 0x11u,    /* set failed (invalid parameter)      */
  UNUR_ERR_DISTR_GET      = 0x12u,    /* get failed (parameter not set)      */
  UNUR_ERR_DISTR_NPARAMS  = 0x13u,    /* invalid number of parameters        */
  UNUR_ERR_DISTR_DOMAIN   = 0x14u,    /* parameter out of domain             */
  UNUR_ERR_DISTR_GEN      = 0x15u,    /* invalid variant for special generator */
  UNUR_ERR_DISTR_REQUIRED = 0x16u,    /* incomplete distribution object, entry missing */
  UNUR_ERR_DISTR_UNKNOWN  = 0x17u,    /* unknown distribution, cannot handle */
  UNUR_ERR_DISTR_INVALID  = 0x18u,    /* invalid distribution object         */

  /** parameter object **/
  UNUR_ERR_PAR_SET        = 0x21u,    /* set failed (invalid parameter)      */
  UNUR_ERR_PAR_VARIANT    = 0x22u,    /* invalid variant -> using default    */
  UNUR_ERR_PAR_INVALID    = 0x23u,    /* invalid parameter object            */

  /** generator object **/
  UNUR_ERR_GEN            = 0x31u,    /* bit for generator object            */
  UNUR_ERR_GEN_DATA       = 0x32u,    /* (possible) invalid data             */
  UNUR_ERR_GEN_CONDITION  = 0x33u,    /* condition for method violated       */
  UNUR_ERR_GEN_INVALID    = 0x34u,    /* invalid generator object            */

  /** misc **/
  UNUR_ERR_ROUNDOFF       = 0x01u,    /* (serious) round-off error           */
  UNUR_ERR_MALLOC         = 0x02u,    /* virtual memory exhausted            */
  UNUR_ERR_NULL           = 0x03u,    /* invalid NULL pointer                */ 
  UNUR_ERR_COOKIE         = 0x04u,    /* invalid cookie                      */
  UNUR_ERR_GENERIC        = 0x05u,    /* generic error                       */

  /** compilation switches **/
  UNUR_ERR_COMPILE        = 0x0eu,    /* not available, recompile library    */

  /** this should not happen **/
  UNUR_ERR_SHOULD_NOT_HAPPEN = 0x0fu, /* error should not happen, report this! */

};

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_ERRNO_H_SEEN */
/*---------------------------------------------------------------------------*/
