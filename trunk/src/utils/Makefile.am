## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/src

noinst_LIBRARIES = libutils.a

libutils_a_SOURCES = \
 debug.c  debug.h debug_source.h \
 error.c  \
 umath.c  umath.h umath_source.h \
 slist.c  slist.h slist_struct.h \
 stream.c         stream_source.h \
 umalloc.c        umalloc_source.h

noinst_HEADERS = \
 unur_fp_source.h \
 unur_errno.h \
 unur_error_source.h \
 unur_math_source.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

