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
#ifndef CODEGEN_SOURCE_H_SEEN
#define CODEGEN_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_stddistr.h>
#include <parser/functparser_source.h>
#include <unuran_acg.h>
#include <time.h>

/*---------------------------------------------------------------------------*/

int _unur_acg_C_header       ( FILE *out, const UNUR_DISTR *distr, const char *rand );
int _unur_acg_FORTRAN_header ( FILE *out, const UNUR_DISTR *distr, const char *rand );
int _unur_acg_JAVA_header    ( FILE *out, const UNUR_DISTR *distr, const char *rand );
int _unur_acg_UNURAN_header  ( FILE *out, const UNUR_DISTR *distr, const char *rand );
/*---------------------------------------------------------------------------*/
/* Code generator for file header.                                           */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_PDF       ( FILE *out, const UNUR_DISTR *distr, const char *pdf );
int _unur_acg_FORTRAN_PDF ( FILE *out, const UNUR_DISTR *distr, const char *pdf );
int _unur_acg_JAVA_PDF    ( FILE *out, const UNUR_DISTR *distr, const char *pdf );
int _unur_acg_UNURAN_PDF  ( FILE *out, const UNUR_DISTR *distr, const char *pdf );
int _unur_acg_UNURAN_PDFbody ( FILE *out, const UNUR_DISTR *distr );
/*---------------------------------------------------------------------------*/
/* Code generator for PDFs of UNURAN build-in standard distributions.        */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_tdr_ps       ( FILE *out, const UNUR_GEN *gen,
			       const char *rand_name, const char *pdf_name );
int _unur_acg_FORTRAN_tdr_ps ( FILE *out, const UNUR_GEN *gen,
			       const char *rand_name, const char *pdf_name );
int _unur_acg_JAVA_tdr_ps    ( FILE *out, const UNUR_GEN *gen,
			       const char *rand_name, const char *pdf_name );
int _unur_acg_UNURAN_tdr_ps  ( FILE *out, const UNUR_GEN *gen,
			       const char *rand_name, const char *pdf_name,
			       int n_cpoints );
/*---------------------------------------------------------------------------*/
/* Code generator for method TDR variant PS (proportional squeeze).          */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_demo_urng       ( FILE *out );
int _unur_acg_FORTRAN_demo_urng ( FILE *out );
int _unur_acg_JAVA_urng         ( FILE *out );
/*---------------------------------------------------------------------------*/
/* Uniform random number generator (for demo mode only).                     */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_main       ( FILE *out, const char *rand_name );
int _unur_acg_FORTRAN_main ( FILE *out, const char *rand_name );
int _unur_acg_JAVA_main    ( FILE *out, const char *class_name );
/*---------------------------------------------------------------------------*/
/* Print main                                                                */
/*---------------------------------------------------------------------------*/

int _unur_acg_C_print_section_rule  ( FILE *out );
int _unur_acg_C_print_section_line  ( FILE *out, const char *format, ... );
int _unur_acg_C_print_section_title ( FILE *out, const char *title );

int _unur_acg_FORTRAN_print_section_rule  ( FILE *out );
int _unur_acg_FORTRAN_print_section_line  ( FILE *out, const char *format, ... );
int _unur_acg_FORTRAN_print_section_title ( FILE *out, const char *title );

int _unur_acg_JAVA_print_section_rule  ( FILE *out );
int _unur_acg_JAVA_print_section_line  ( FILE *out, const char *format, ... );
int _unur_acg_JAVA_print_section_title ( FILE *out, const char *title );
/*---------------------------------------------------------------------------*/
/* Print a section header with n_lines lines to output stream.               */
/*---------------------------------------------------------------------------*/

int _unur_acg_FORTRAN_print_double ( FILE *out, double x );
/*---------------------------------------------------------------------------*/
/* Print a double number in FORTRAN format.                                  */
/*---------------------------------------------------------------------------*/

int _unur_acg_JAVA_begin_class ( FILE *out, const UNUR_GEN *gen, const char *class ); 
int _unur_acg_JAVA_end_class   ( FILE *out, const UNUR_GEN *gen );
/*---------------------------------------------------------------------------*/
/* JAVA: begin and end of class containing PDF and tdr                       */
/*---------------------------------------------------------------------------*/

int _unur_acg_JAVA_tdr_class_IV( FILE *out, const UNUR_GEN *gen );
/*---------------------------------------------------------------------------*/
/* JAVA: definition of class IV                                              */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif   /* CODEGEN_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
