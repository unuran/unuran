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

/*

=NODE  StringAPI   String Interface
=UP TOP [25]

=DESCRIPTION

   The string interface provided by the @command{unur_str2gen} call is
   the easiest way to use UNURAN. The function takes a character
   string as its argument. The string is parsed and the information
   obtained is used to create a generator object. It returns NULL if
   this fails, either due to a syntax error, or due to invalid data.
   Notice that the string interface does not implement all features of
   the UNURAN library. For trickier tasks it might be necessary to use
   the UNURAN API. Especially using generic distributions is not fully
   supported yet.

=END

*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_GEN *unur_str2gen( const char *string );
/* 
   Get a generator object for the distribution, method and unifrom
   random number generator as discribed in the given @var{string}.
   See @ref{StringSyntax,Syntax of String Interface,Syntax of String
   Interface} for details.
*/

/*
=EON
*/

/*---------------------------------------------------------------------------*/
/*

=NODE  StringSyntax   Syntax of String Interface
=UP StringAPI [10]

=DESCRIPTION

   The given string holds information about the requested distribution
   and (optional) about the sampling method and the uniform random
   number generator invoked. The interpretation of the string is not
   case-sensitive, all white spaces are ignored.

   The string consists of up to three blocks, separated by colons
   @code{:}.

   Each block consists of @code{<key>=<value>} pairs, separated by
   semicolons @code{;}.

   The first key in each block is used to indicate each block.

   We have three different blocks with the following (first) keys:
   @table @code 
   @item distr
      definition of the distribution (see 
      @ref{StringDistr, Distribution String, Distribution String}).

   @item method
      description of the transformation method
      (see @ref{StringMethod, Method String, Method String}).

   @item urng
      uniform random number generation
      (see @ref{StringURNG,Uniform RNG String,Uniform RNG String}).
   @end table

   The @code{distr} block must be the very first block and is
   obligatory. All the other blocks are optional and can be arranged
   in arbitrary order.

   For details see the following description of each block.

   The other @code{<key>=<value>} pairs are used to set parameters.
   The name of the parameter is given by the @code{<key>} string. It is
   deduced from the UNURAN set calls by taking the part after
   @code{@dots{}_set_}.
   The @code{<value>} string holds tokens for the parameters to be
   set, separated by commata @code{,}.
   There are two types of tokens:
   @table @emph
   @item single token
      that represents a number, and

   @item list
      i.e., a list of single tokens, separated by commata @code{,},
      enclosed in parenthesis @code{(...)}.
   @end table

   The @code{<value>} string (including the character @code{=}) can be
   omitted when no argument is required.

   At the moment not all @command{set} calls are supported.
   The syntax for the @code{<value>} can be directly derived from the
   corresponding @command{set} calls. To simplify the syntax additional
   shortcuts are possible. The following table list the parameters for
   the @code{set} calls that are supported by the string interface; the
   entry in parenthesis gives the type of the argument in the string:

   @table @code
   @item int  @r{(}token@r{):}
      The token is interpreted as an integer number.
      @code{true} and @code{on} are transformed to @code{1},
      @code{false} and @code{off} are transformed to @code{0}.
      A missing argument is interpreted as @code{1}.

   @item unsigned @r{(}token@r{):} 
      The token is interpreted as an unsigned hexadecimal integer number.

   @item double  @r{(}token@r{):}
      The token is interpreted as a floating point number.

   @item double, double  @r{(}token, token @r{or} list@r{):}
      The two tokens or the first two entries in the list are interpreted as
      a floating point numbers.

   @item int, double*  @r{(}[token,] list @r{or} token@r{):}
      @itemize @minus
      @item
         The list is interpreted as a double array.
	 The (first) token as its length.
	 If it is less than the actual size if the array only the first entries
	 of the array is used.
      @item
         If only the list is given (i.e., if the first token is omitted),
	 it is replaced by the actual size of the array.
      @item
	 If only the token is given (i.e., if the list is omitted), the NULL
	 pointer is used instead an array as argument.
      @end itemize

   @item double*, int  @r{(}list [,token]@r{):}
      The list is interpreted as a double array.
      The (second) token as its length. 
      If the length is omitted, it is replaced by the actual size of the
      array. (Only in the @code{distribution} block!)

   @end table

=EON */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*

=NODE  StringDistr    Distribution String
=UP StringAPI [20]

=DESCRIPTION

   The @code{distr} block must be the very first block and is
   obligatory. For that reason the keyword @code{distr} is optional and
   can be omitted. Moreover it is ignored while parsing the string. To
   avoid some possible confusion, however, it has to start with the
   letter @code{d} (if it is given at all). 

   The value of the @code{distr} key is used to get the distribution
   object, either via a @command{unur_distr_<value>} call for a standard
   distribution. The parameters for the standard distribution are given
   as a list. There must not be any character (other than white space)
   between the name of the standard distribution and the opening
   parenthesis @code{(} of this list. E.g., to get a beta distribution,
   use
   @smallexample
      distr = beta(2,4)
   @end smallexample

   Or via a @command{unur_distr_<value>_new} call to get an object of a generic
   distribution. However not all generic distributions are supported yet.
   E.g., to get an object for a discrete distribution with probability
   vector (1,2,3), use 
   @smallexample
      distr = discr; pv = (1,2,3)
   @end smallexample

=EON
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*

=NODE  StringMethod   Method String
=UP StringAPI [30]

=DESCRIPTION

   The key @code{method} is obligatory, it must be the first key and its
   value is the name of a method suitable for the choosen standard
   distribution. E.g., if method AROU is chosen, use
   @smallexample
      method = arou
   @end smallexample

   Of course the all following keys dependend on the method choosen at
   first. All corresponding @command{set} calls of UNURAN are available
   and the key is the string after the @command{unur_<methodname>_set_}
   part of the command. E.g., UNURAN provides the command 
   @command{unur_arou_set_max_sqhratio} to set a parameter of method AROU.
   To call this function via the string-interface, the
   key @code{max_sqhratio} can be used:
   @smallexample
      max_sqhratio = 0.9
   @end smallexample
   
   If this block is omitted, a suitable default method is used. Notice
   however that the default method may change in future versions of
   UNURAN. Moreover, only standard distributions are supported yet.

=EON
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*

=NODE  StringURNG    Uniform RNG String
=UP StringAPI [40]

=DESCRIPTION

   The value of the @code{urng} key is passed to the PRNG interface (see
   @ifinfo
      @xref{Top, , Overview, prng, PRNG Manual}.
   @end ifinfo
   @ifnotinfo
      @uref{http://statistik.wu-wien.ac.at/prng/manual/,PRNG manual}
   @end ifnotinfo
   for details).
   However it only works when using the PRNG library is enabled, 
   see @ref{Installation} for details. There are no other keys.

   If this block is omitted the UNURAN default generator is used.

=EON
*/
/*---------------------------------------------------------------------------*/

/*
Example:

The standard distribution block must be first!
Here an example to illustrate that, the details will be explained below:
unur_str2gen("distr=normal(0,1) : method=arou; max_sqhratio=0.9 : prng=MT19937(133)")
*/

/*---------------------------------------------------------------------------*/





