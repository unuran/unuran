/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_set.h                                                        *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros for indicating which parameters have been set.     *
 *         defines function prototypes for setting, changing or reading      *
 *         parameters in generator objects.                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only used in unur_methods.h                                       *
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
#ifndef __UNUR_SET_H_SEEN
#define __UNUR_SET_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*  Mark parameter settings                                                  */

#define UNUR_SET_DOMAIN      0x00000001UL
#define UNUR_SET_PDFPARAM    0x00000002UL
#define UNUR_SET_MODE        0x00000004UL
#define UNUR_SET_AREA        0x00000008UL

#define UNUR_SET_CENTER      0xf0000001UL

#define UNUR_SET_VARIANT     0x00000010UL
#define UNUR_SET_FACTOR      0x00000020UL

#define UNUR_SET_STP         0x00000100UL
#define UNUR_SET_N_STP       0x00000200UL
#define UNUR_SET_MAX_RATIO   0x00000400UL
#define UNUR_SET_MAX_IVS     0x00000800UL

#define UNUR_SET_TDR_C       0x10001000UL
#define UNUR_SET_TABL_C      0x10002000UL
#define UNUR_SET_SLOPES      0x10004000UL

/*---------------------------------------------------------------------------*/
/* Parameters for the distribution and its p.d.f.                            */

int unur_set_domain( struct unur_par *parameter, double left, double right );
/* set the left and right borders of the domain of the distribution          */

int unur_set_domain_vec( struct unur_par *par, double **domain );
/* set coordinates for domain boundary                                       */

int unur_set_mode( struct unur_par *parameter, double mode );
/* set mode of p.d.f.                                                        */

int unur_set_usemode( struct unur_par *par, int usemode );
/* set flag for using mode of p.d.f.                                         */

/*---------------------------------------------------------------------------*/
/* Parameters for generators of univariate discrete distributions            */

int unur_set_factor( struct unur_par *parameter, double factor );
/* set factor for relative size of (search|guide|alias) table                */


/*---------------------------------------------------------------------------*/
/* Parameters for generators of univariate continuous distributions          */
 
int unur_set_cpoints( struct unur_par *parameter, int n_stp, double *starting_cpoints );
/* set construction points for hat and/or its number for initialization      */

int unur_get_n_intervals( struct unur_gen *gen );
/* get number of intervals/segments                                          */

int unur_get_max_intervals( struct unur_gen *gen );
/* get maximal number of intervals/segments                                  */

int unur_set_max_shratio( struct unur_par *parameter, double max_ratio );
/* set bound for ratio A(squeeze) / A(hat)                                   */

double unur_get_shratio( struct unur_gen *gen );
/* get ratio A(squeeze) / A(hat)                                             */

int unur_set_max_intervals( struct unur_par *parameter, int max_ivs );
/* set maximum number of intervals or segments                               */

int unur_set_tdr_c( struct unur_par *par, double c );
/* set parameter c for transformation T_c  (method TDR only)                 */           

int unur_set_tabl_c( struct unur_par *par, double c );
/* set parameter for equal area rule  (method TABL only)                     */           

int unur_set_slopes( struct unur_par *par, double *slopes, int n_slopes );
/* set slopes of p.d.f. (method TABL only)                                   */


/*---------------------------------------------------------------------------*/
/* Parameters for generators of mulitvariate continuous distributions        */
int unur_get_dimension( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* Parameters for all generators                                             */

int unur_set_urng( struct unur_par *par, UNUR_URNG_TYPE urng );
/* set uniform random number generator                                       */

UNUR_URNG_TYPE unur_chg_urng( struct unur_gen *gen, UNUR_URNG_TYPE urng );
/* change uniform random number generator                                    */

UNUR_URNG_TYPE unur_get_urng( struct unur_gen *gen );
/* get uniform random number generator                                    */

int unur_set_variant( struct unur_par *par, unsigned long variant );
/* set variant of method                                                     */

int unur_set_check(  struct unur_par *parameter, int check );
/* turn testing of sampling on/off */

int unur_set_copyall(  struct unur_par *par, int copy );
/* turn copaing of all inputs into generator object on/off                   */

int unur_set_debug( struct unur_par *parameter, unsigned long db );
/* set debugging flag for generator                                          */

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_SET_H_SEEN */
/*---------------------------------------------------------------------------*/

