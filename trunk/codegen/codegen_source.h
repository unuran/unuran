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
#include "PDFgen_source.h"

/*---------------------------------------------------------------------------*/

int _unur_tdr_ps_codegen( struct unur_gen *gen, FILE *out, 
			  const char *rand_name, const char *pdf_name );
/*---------------------------------------------------------------------------*/
/* Code generator for method TDR variant PS (proportional squeeze).          */
/*---------------------------------------------------------------------------*/

void _unur_acg_print_sectionheader( FILE *out, int n_lines, ... );
/*---------------------------------------------------------------------------*/
/* Print a section header with n_lines lines to output stream.               */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif   /* __CODEGEN_H_SEEN */
/*---------------------------------------------------------------------------*/
