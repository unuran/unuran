## Process this file with automake to produce Makefile.in
## $Id$

INCLUDES = -I$(top_srcdir)/src

noinst_LTLIBRARIES = libutils.la

libutils_la_SOURCES = \
 debug.c  debug.h debug_source.h \
 eigensystem.c \
 error.c  \
 fmax.c 	  fmax_source.h \
 hooke.c 	  hooke_source.h \
 matrix.c 	  matrix_source.h \
 slist.c  slist.h slist_struct.h \
 stream.c         stream_source.h \
 string.c         string_source.h string_struct.h \
 umalloc.c        umalloc_source.h \
 umath.c  umath.h umath_source.h \
 unur_fp.c

noinst_HEADERS = \
 unur_fp_source.h unur_fp_const_source.h \
 unur_errno.h \
 unur_error_source.h \
 unur_math_source.h

# clean backup files
CLEANFILES = *~

# clean generated files
##DISTCLEANFILES = \
## unur_fp_const_source.h

MAINTAINERCLEANFILES = \
 Makefile.in \
 unur_fp_const_source.h

