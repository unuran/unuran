/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: functparser_struct.h                                              *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for function parser                           *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
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
/* ?????? */
#define SYMBLENGTH 20 

/*---------------------------------------------------------------------------*/
/* Structure for function tree                                               */

struct treenode { 
  char            symb[SYMBLENGTH];  /* zeigt auf Symbol aus Symboltab. */ 
  int             token;             /* Token des Symbols               */ 
  int             symbkind;          /* Art des Symbols (REL_OP etc.)   */ 
  float           val;               /* aktueller arithmetischer Wert   */ 
  struct treenode *left;             /* Zeiger auf linken Sohn          */ 
  struct treenode *right;            /* Zeiger auf rechten Sohn         */ 

#ifdef UNUR_COOKIES
  unsigned cookie;              /* magic cookie                              */
#endif
}; 

/*---------------------------------------------------------------------------*/
/* Symbols used in function string                                           */

struct symbols { 
  char            name[SYMBLENGTH];  /* Name des Symbols (z. B. "SIN")  */ 
  int             info;              /* Prioritaet bzw. Argumentanzahl  */ 
  double           val;               /* Konstanten: numerischer Wert    */ 
  double           (*vcalc)(int t, double l, double r);        
                                     /* Zeiger auf Berechnungsfunktion  */ 
  char            *(*dcalc)(char *par,struct treenode *w,
                            char *l, char *r, char *dl, char *dr,char *s);
                                     /* Zeiger auf Ableitungsfunktion   */ 
  struct treenode *tree;             /* Bei UFUNCS: Zeiger auf Baum     */ 

};

/*---------------------------------------------------------------------------*/
