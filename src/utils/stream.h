/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: stream.h                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         routines for output streams and reading data                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
   =NODE  Output_streams Output streams

   =UP Errno

   =DESCRIPTION
      @cindex Error handlers
      @cindex Output streams

      In addition to reporting error via the @code{unur_errno}
      mechanism UNU.RAN writes a short error diagnostics to an output
      stream, usually an open file handler. This stream can be set at
      runtime by the unur_set_stream() call. If no such stream is
      given by the user a default stream is used by the library: all
      warnings and error messages are written into the file
      @file{unuran.log} in the current working directory. The name of
      this log file is defined by the macro @code{UNUR_LOG_FILE} in 
      @file{unuran_config.h}. 
      
      This output stream is also used to log descriptions of built
      generator objects and for writing debugging information.
      If you want to use this output stream for your own programs use 
      unur_get_stream() to get its file handler.
      
      All warnings, error messages and all debugging information
      are written onto the same output stream.
      To destinguish between the messages for different generators 
      every generator object has its own identifier which is 
      composed by the generator type, followed by a dot and three digits.
      (If there are more than 999 generators then the identifiers are
      not unique.)

   =END      
*/

/* =ROUTINES */

/*---------------------------------------------------------------------------*/
/* manipulate output stream                                                  */

FILE *unur_set_stream( FILE *new_stream );
/*
  Set new file handle for output stream; the old file handle is
  returned. The NULL pointer is not allowed. (If you want to disable
  logging of debugging information use 
  unur_set_default_debug(UNUR_DEBUG_OFF) instead.)

  The output stream is used to report errors and warning, and
  debugging information. It is also used to log descriptions of
  build generator objects (when this feature is switched on; see also ?).
*/

FILE *unur_get_stream( void );
/*
  Get the file handle for the current output stream.
*/

/* =END */

/*---------------------------------------------------------------------------*/
