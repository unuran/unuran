/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: parser.h                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for parser                                    *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_GEN *unur_str2gen( const char *string );
/* 
   Get ....
*/





/********************************************************************/
/*

 *   There are additional transformations before executing the necessary     *
 *   set calls (listed by argument types)                                    *
 *      int:      'true' and 'on' are transformed to 1,                      *
 *                'false' and 'off' are transformed to 0.                    *
 *                a missing argument is transformed to 1.                    *
 *      unsigned: string is interpreted as hexadecimal number.               *
 *      double:   none                                                       *
 *      double, double: Instead of two doubles, a list with at least two     *
 *                numbers can be given.                                      *
 *      int, double*: If the first argument is missing, it is replaced       *
 *                by the size of the array.                                  *
 *                If the second argument is missing, the NULL pointer is     *
 *                used instead an array as argument.                         *
 *      double*, int: Only for distributions!                                *
 *                If the second argument 'int' is missing, it is replaced    *
 *                by the size of the array.                                  *

 */
/********************************************************************/
/*

String Interface for .....


The string interface provided by the unur_str2gen() call is the
easiest way to use UNURAN. The function takes a character string as
its argument. The string is parsed and the information obtained is
used to create a generator object. It return NULL if this fails,
either due to a syntax error, or due to invalid data.
Notice that the string interface does not implement all features of
the UNURAN library. For trickier tasks it might be necessary to use
the UNURAN API. Especially using generic distributions is not fully 
supported yet.


The string interface:
The given string holds information about the requested distribution
and (optional) about the sampling method and the
uniform random number generator invoked.
The syntax of the string is not case-sensitive, all white spaces are
ignored.
The string consists of up to three blocks, separated by colons
@code{:}. Each block consists of <key>=<value> pairs, separated by
semicolons @code{;}. 
The first key in each block is used to indicate
each block.
The other <key>=<value> pairs are used to set parameters.
The name of the parameter is given by the <key> string.
It is deduced from the UNURAN set calls by taking the part after
"..._set_".
The <value> string holds tokens for the parameters to be set,
separated by commata @code{,}.
There are two types of tokens:
      single tokens, that represent numbers, and
      list, i.e. a list of single tokens, separated by commata @code{,},
            enclosed in parenthesis @code{(...)}.
The <value> string (including the character '=') can be omitted when
no argument is required.


We have three different blocks with keys
    distr  ... contains the definition of the distribution
    method ... contains the description of the transformation method
    urng   ... contains the definition of the uniform RNG

The 'distr' block must be the very first block and is obligatory.
All the other blocks are optional and can be arranged in arbitrary order.

Distribution:
The 'distr' block must be the very first block and is obligatory.
For that reason the keyword 'distr' is optional and can be omitted.
Moreover it is ignored while parsing the string. To avoid some
possible confussion it only has to start with the letter @code{d} 
(if it is given at all.) 
The value of the 'distr' key is used to get the distribution object,
either via a unur_distr_<value>() call for a standard distribution.
The parameters for the standard distribution are given as a list.
There must not be any character (other than white space) between the
name of the standard distribution and the opening parenthesis @code{(}
of this list.
Or via a unur_distr_<value>_new() call to get an object of a generic
distribution. However not all generic distributions are supported yet.
E.g.: distr = beta(2,4)



Method:
The key `method' is obligatory, it must be first and the value is the
name of a method suitable for the choosen standard distribution.
E.g.: method = arou
If this block is omitted, a suitable default method is used. Notice
however that the default method may change in following versions of
UNURAN.
Of course the following keys dependend on the method choosen at first.
All corresponding `_set_' commands of UNURAN are available and the key
is the string after the `unur_<methodname>_set_' part of the command.
e.g.: UNURAN provides the command `unur_arou_set_max_sqhratio' to
set a parameter of the method AROU.
To call this function via the string-interface, the
key `max_sqhratio' can be used:
max_sqhratio = 0.9.


Uniform random number generator:
The value of the 'urng' key is passed to the PRNG interface (see the
PRNG manual for details).
However it only works when using the PRNG library is enabled, i.e.
the libray must be installed and the corresponding
compiler switch must be set in @file{unuran_config.h}.
If this block is omitted the UNURAN default generator is used.


Example:



The standard distribution block must be first!
Here an example to illustrate that, the details will be explained below:
unur_str2gen("distr=normal(0,1) : method=arou; max_sqhratio=0.9 : prng=MT19937(133)")


 */


/* =END */
/*---------------------------------------------------------------------------*/

