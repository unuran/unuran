/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      codegen_source.h                                             *
 *                                                                           *
 *   internal declarations for code generators                               *
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
#ifndef __CODEGEN_H_SEEN
#define __CODEGEN_H_SEEN
/*---------------------------------------------------------------------------*/

#include <source_unuran.h>
#include <time.h>

/*---------------------------------------------------------------------------*/

int _unur_acg_C_header (UNUR_DISTR *distr, FILE *out, const char *rand);
int _unur_acg_FORTRAN_header (UNUR_DISTR *distr, FILE *out, const char *rand);
int _unur_acg_JAVA_header (UNUR_DISTR *distr, FILE *out, const char *rand);
int _unur_acg_UNURAN_header (UNUR_DISTR *distr, FILE *out, const char *rand);
/*---------------------------------------------------------------------------*/
/* Code generator for file header.                                           */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf);
int _unur_acg_FORTRAN_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf);
int _unur_acg_JAVA_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf);
/*---------------------------------------------------------------------------*/
/* Code generator for PDFs of UNURAN build-in standard distributions.        */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_tdr_ps( struct unur_gen *gen, FILE *out, 
			const char *rand_name, const char *pdf_name );
int _unur_acg_FORTRAN_tdr_ps( struct unur_gen *gen, FILE *out, 
			      const char *rand_name, const char *pdf_name );
int _unur_acg_JAVA_tdr_ps( struct unur_gen *gen, FILE *out, 
			   const char *rand_name, const char *pdf_name );
/*---------------------------------------------------------------------------*/
/* Code generator for method TDR variant PS (proportional squeeze).          */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_demo_urng( FILE *out );
int _unur_acg_FORTRAN_demo_urng( FILE *out );
int _unur_acg_JAVA_urng( FILE *out );
/*---------------------------------------------------------------------------*/
/* Uniform random number generator (for demo mode only).                     */
/*---------------------------------------------------------------------------*/

void _unur_acg_C_print_section_rule( FILE *out );
void _unur_acg_C_print_section_line( FILE *out, const char *format, ... );
void _unur_acg_C_print_section_title( FILE *out, const char *title );

void _unur_acg_FORTRAN_print_section_rule( FILE *out );
void _unur_acg_FORTRAN_print_section_line( FILE *out, const char *format, ... );
void _unur_acg_FORTRAN_print_section_title( FILE *out, const char *title );

void _unur_acg_JAVA_print_section_rule( FILE *out );
void _unur_acg_JAVA_print_section_line( FILE *out, const char *format, ... );
void _unur_acg_JAVA_print_section_title( FILE *out, const char *title );
/*---------------------------------------------------------------------------*/
/* Print a section header with n_lines lines to output stream.               */
/*---------------------------------------------------------------------------*/

void _unur_acg_FORTRAN_print_double( FILE *out, double x );
/*---------------------------------------------------------------------------*/
/* Print a double number in FORTRAN format.                                  */
/*---------------------------------------------------------------------------*/

int _unur_acg_JAVA_begin_class( struct unur_gen *gen, FILE *out, const char *class ); 
int _unur_acg_JAVA_end_class( struct unur_gen *gen, FILE *out );
/*---------------------------------------------------------------------------*/
/* JAVA: begin and end of class containing PDF and tdr                       */
/*---------------------------------------------------------------------------*/

int _unur_acg_JAVA_tdr_class_IV( struct unur_gen *gen, FILE *out );
/*---------------------------------------------------------------------------*/
/* JAVA: definition of class IV                                              */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif   /* __CODEGEN_H_SEEN */
/*---------------------------------------------------------------------------*/
