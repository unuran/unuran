## Process this file with automake to produce Makefile.in
## $Id$

AM_CFLAGS = -Wall

INCLUDES = \
 -I$(top_srcdir)/src \
 -I$(top_srcdir)/src/utils \
 -I$(top_srcdir)/src/methods

noinst_LIBRARIES = libstdgen.a

libstdgen_a_SOURCES = \
 exponential.c \
 gamma.c \
 normal.c

pkginclude_HEADERS = \
 unur_stdgen.h

noinst_HEADERS = \
 unur_stdgen_variants.h

# clean backup files too
clean-local:
	rm -f *~
